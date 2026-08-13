#!/usr/bin/env python3
"""Build the static Reveal.js thesis deck from its PowerPoint source."""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import posixpath
import shutil
import struct
import subprocess
import tempfile
import zipfile
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from xml.etree import ElementTree as ET


PRESENTATION_NS = "http://schemas.openxmlformats.org/presentationml/2006/main"
DRAWING_NS = "http://schemas.openxmlformats.org/drawingml/2006/main"
RELATIONSHIP_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PACKAGE_RELATIONSHIP_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
NS = {"p": PRESENTATION_NS, "a": DRAWING_NS, "r": RELATIONSHIP_NS}

EXPECTED_WIDTH = 12192000
EXPECTED_HEIGHT = 6858000
EXPECTED_SOURCE_SHA256 = "7bf130f07bc057d275ef2d5a5403edde3b8efa632d5558c0052c9479cf6becc9"
RENDERER_IMAGE = "craftgrn-presentation-renderer:ubuntu24.04-lo24.2.7-poppler24.02.0"
RENDERER_PROVENANCE = {
    "mode": "docker",
    "base_image": "ubuntu@sha256:561618e2c15bf2397621dd04f96926663a3b5616c189cf7e38db7e82f5c538ea",
    "ubuntu_snapshot": "20260813T001500Z",
    "certificate_bootstrap": "noble-release:ca-certificates=20240203,openssl=3.0.13-0ubuntu3",
    "libreoffice": "4:24.2.7-0ubuntu0.24.04.6",
    "poppler_utils": "24.02.0-1ubuntu9.9",
}

DEMOS = {
    14: {
        "url": "https://oncologylab.github.io/fp-tools/demos/reports/diff_footprints_K562_HepG2.html",
        "label": "Activate differential footprint report",
        "left": 4.747,
        "top": 9.494,
        "width": 90.506,
        "height": 90.506,
    },
    20: {
        "url": "https://oncologylab.github.io/fp-tools/demos/gui/fp-tools-gui-static-demo.html#home",
        "label": "Activate fp-tools GUI demo",
        "left": 14.611,
        "top": 15.351,
        "width": 69.894,
        "height": 75.132,
    },
    43: {
        "url": "demos/IRF4_Fibroblast_top100_log2fc0p5_network.html",
        "label": "Activate interactive IRF4 network",
        "left": 23.436,
        "top": 15.450,
        "width": 53.145,
        "height": 75.547,
    },
}


@dataclass(frozen=True)
class Slide:
    source_number: int
    title: str
    notes_text: str
    notes_html: str


@dataclass(frozen=True)
class Presentation:
    width: int
    height: int
    slides: list[Slide]


def _relationship_path(source: str) -> str:
    path = PurePosixPath(source)
    return str(path.parent / "_rels" / f"{path.name}.rels")


def _relationship_map(archive: zipfile.ZipFile, source: str) -> dict[str, tuple[str, str]]:
    relationship_path = _relationship_path(source)
    if relationship_path not in archive.namelist():
        return {}
    root = ET.fromstring(archive.read(relationship_path))
    return {
        node.attrib["Id"]: (
            node.attrib.get("Type", "").rsplit("/", 1)[-1],
            node.attrib["Target"],
        )
        for node in root
    }


def _resolve_relationship(source: str, target: str) -> str:
    return posixpath.normpath(posixpath.join(posixpath.dirname(source), target))


def _ascii_html(text: str) -> str:
    escaped = html.escape(text, quote=True)
    return escaped.encode("ascii", "xmlcharrefreplace").decode("ascii")


def _shape_text(shape: ET.Element) -> str:
    paragraphs = []
    for paragraph in shape.findall(".//a:p", NS):
        value = "".join(node.text or "" for node in paragraph.findall(".//a:t", NS))
        if value:
            paragraphs.append(value)
    return "\n".join(paragraphs)


def _slide_title(root: ET.Element) -> str:
    shapes = root.findall("./p:cSld/p:spTree/p:sp", NS)
    for shape in shapes:
        placeholder = shape.find("./p:nvSpPr/p:nvPr/p:ph", NS)
        if placeholder is not None and placeholder.attrib.get("type") in {"title", "ctrTitle"}:
            return _shape_text(shape).replace("\n", " ").strip()
    for shape in shapes:
        value = _shape_text(shape).replace("\n", " ").strip()
        if value:
            return value
    return ""


def _run_html(run: ET.Element) -> tuple[str, str]:
    text_node = run.find("a:t", NS)
    text = text_node.text if text_node is not None and text_node.text is not None else ""
    value = _ascii_html(text)
    properties = run.find("a:rPr", NS)
    if properties is not None:
        if properties.attrib.get("b") in {"1", "true"}:
            value = f"<strong>{value}</strong>"
        if properties.attrib.get("i") in {"1", "true"}:
            value = f"<em>{value}</em>"
    return text, value


def _notes_content(root: ET.Element) -> tuple[str, str]:
    text_paragraphs = []
    html_paragraphs = []
    for shape in root.findall("./p:cSld/p:spTree/p:sp", NS):
        placeholder = shape.find("./p:nvSpPr/p:nvPr/p:ph", NS)
        placeholder_type = placeholder.attrib.get("type", "body") if placeholder is not None else "body"
        if placeholder_type in {"sldImg", "sldNum", "hdr", "ftr", "dt"}:
            continue
        for paragraph in shape.findall("./p:txBody/a:p", NS):
            text_parts = []
            html_parts = []
            for child in paragraph:
                local_name = child.tag.rsplit("}", 1)[-1]
                if local_name in {"r", "fld"}:
                    plain, marked_up = _run_html(child)
                    text_parts.append(plain)
                    html_parts.append(marked_up)
                elif local_name == "br":
                    text_parts.append("\n")
                    html_parts.append("<br>")
            plain_paragraph = "".join(text_parts)
            marked_up_paragraph = "".join(html_parts)
            if plain_paragraph or marked_up_paragraph:
                text_paragraphs.append(plain_paragraph)
                html_paragraphs.append(f"<p>{marked_up_paragraph}</p>")
    return "\n".join(text_paragraphs), "".join(html_paragraphs)


def parse_presentation(source: Path) -> Presentation:
    with zipfile.ZipFile(source) as archive:
        presentation_path = "ppt/presentation.xml"
        root = ET.fromstring(archive.read(presentation_path))
        size = root.find("p:sldSz", NS)
        if size is None:
            raise ValueError("PowerPoint source does not define a slide size")
        width = int(size.attrib["cx"])
        height = int(size.attrib["cy"])
        if width * 9 != height * 16:
            raise ValueError(f"PowerPoint source must use a 16:9 slide size; found {width}x{height} EMU")

        presentation_relationships = _relationship_map(archive, presentation_path)
        slides = []
        for slide_id in root.findall("p:sldIdLst/p:sldId", NS):
            relationship_id = slide_id.attrib[f"{{{RELATIONSHIP_NS}}}id"]
            _, target = presentation_relationships[relationship_id]
            slide_path = _resolve_relationship(presentation_path, target)
            slide_root = ET.fromstring(archive.read(slide_path))
            source_number = int(PurePosixPath(slide_path).stem.removeprefix("slide"))
            slide_relationships = _relationship_map(archive, slide_path)
            notes_path = None
            for relationship_type, relationship_target in slide_relationships.values():
                if relationship_type == "notesSlide":
                    notes_path = _resolve_relationship(slide_path, relationship_target)
                    break
            if notes_path is None:
                notes_text, notes_html = "", ""
            else:
                notes_root = ET.fromstring(archive.read(notes_path))
                notes_text, notes_html = _notes_content(notes_root)
            slides.append(
                Slide(
                    source_number=source_number,
                    title=_slide_title(slide_root),
                    notes_text=notes_text,
                    notes_html=notes_html,
                )
            )
    return Presentation(width=width, height=height, slides=slides)


def _demo_markup(slide_number: int) -> str:
    demo = DEMOS.get(slide_number)
    if demo is None:
        return ""
    style = (
        f"--demo-left:{demo['left']}%;--demo-top:{demo['top']}%;"
        f"--demo-width:{demo['width']}%;--demo-height:{demo['height']}%"
    )
    url = _ascii_html(str(demo["url"]))
    label = _ascii_html(str(demo["label"]))
    return (
        f'<div class="demo" style="{style}" data-demo-url="{url}">'
        f'<button class="demo-activate" type="button">{label}</button>'
        f'<a class="demo-fallback" href="{url}" target="_blank" rel="noopener">Open demo in new tab</a>'
        '<div class="demo-frame" hidden>'
        '<div class="demo-toolbar">'
        '<span class="demo-status" aria-live="polite">Loading interactive demo...</span>'
        f'<a class="demo-new-tab" href="{url}" target="_blank" rel="noopener">Open in new tab</a>'
        '<button class="demo-exit" type="button">Exit demo</button>'
        '</div><iframe title="Interactive demonstration" loading="lazy" '
        'sandbox="allow-scripts allow-forms allow-downloads allow-popups '
        'allow-popups-to-escape-sandbox"></iframe></div></div>'
    )


def _index_html(presentation: Presentation) -> str:
    sections = []
    for display_number, slide in enumerate(presentation.slides, start=1):
        title = _ascii_html(slide.title or f"Slide {display_number}")
        notes = slide.notes_html
        sections.append(
            f'<section id="slide-{display_number}" data-slide-number="{display_number}" '
            f'data-background-image="assets/slides/slide-{display_number:03d}.png" '
            f'data-background-size="contain" aria-label="{title}">'
            f'{_demo_markup(display_number)}<aside class="notes">{notes}</aside></section>'
        )
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="description" content="Yaoxiang Li thesis presentation">
<title>From chromatin footprints to gene regulatory networks</title>
<link rel="stylesheet" href="vendor/reveal/reveal.css">
<link rel="stylesheet" href="assets/presentation.css">
</head>
<body>
<div class="reveal"><div class="slides">{''.join(sections)}</div></div>
<script src="vendor/reveal/reveal.js"></script>
<script src="vendor/reveal/plugin/notes/notes.js"></script>
<script src="assets/presentation.js"></script>
<script>
Reveal.initialize({{
  width: 1920,
  height: 1080,
  margin: 0,
  minScale: 0.1,
  maxScale: 2,
  controls: true,
  progress: true,
  hash: true,
  history: true,
  center: false,
  transition: "none",
  backgroundTransition: "none",
  slideNumber: "c/t",
  showNotes: false,
  plugins: [ RevealNotes ]
}});
</script>
</body>
</html>
"""


def generate_site(
    presentation: Presentation,
    source: Path,
    output: Path,
    expected_slide_count: int = 59,
    renderer_provenance: dict[str, str] | None = None,
) -> None:
    if len(presentation.slides) != expected_slide_count:
        raise ValueError(
            f"Expected {expected_slide_count} slides, found {len(presentation.slides)}"
        )
    output.mkdir(parents=True, exist_ok=True)
    source_sha256 = hashlib.sha256(source.read_bytes()).hexdigest()
    manifest = {
        "source_filename": source.name,
        "source_sha256": source_sha256,
        "slide_count": len(presentation.slides),
        "slide_size_emu": [presentation.width, presentation.height],
        "render_size_px": [3840, 2160],
        "renderer": renderer_provenance or {"mode": "unspecified"},
        "slides": [
            {
                "number": number,
                "source_number": slide.source_number,
                "title": slide.title,
                "notes_text": slide.notes_text,
                "image_sha256": hashlib.sha256(
                    (output / "assets" / "slides" / f"slide-{number:03d}.png").read_bytes()
                ).hexdigest(),
            }
            for number, slide in enumerate(presentation.slides, start=1)
        ],
    }
    (output / "index.html").write_text(_index_html(presentation), encoding="ascii")
    (output / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=True, indent=2) + "\n",
        encoding="ascii",
    )


RUNTIME_CSS = """html,body{margin:0;width:100%;height:100%;background:#000}
.reveal{background:#000}
.reveal .slides section{box-sizing:border-box;width:100%;height:100%;overflow:hidden}
.reveal .slide-background-content{background-color:#fff}
.demo{position:absolute;left:var(--demo-left);top:var(--demo-top);width:var(--demo-width);height:var(--demo-height);z-index:20}
.demo-activate{position:absolute;right:18px;bottom:18px;padding:13px 18px;border:2px solid #fff;border-radius:6px;background:#17365d;color:#fff;font:700 18px Arial,sans-serif;box-shadow:0 2px 10px rgba(0,0,0,.35);cursor:pointer}
.demo-fallback{position:absolute;right:18px;bottom:72px;padding:8px 12px;border:2px solid #17365d;border-radius:5px;background:#fff;color:#17365d;font:700 15px Arial,sans-serif;text-decoration:none;box-shadow:0 2px 10px rgba(0,0,0,.2)}
.demo-activate:focus-visible,.demo-exit:focus-visible,.demo-new-tab:focus-visible,.demo-fallback:focus-visible{outline:4px solid #f6c344;outline-offset:2px}
.demo-frame{position:absolute;inset:0;background:#fff;border:2px solid #17365d;box-shadow:0 2px 14px rgba(0,0,0,.35)}
.demo-toolbar{box-sizing:border-box;height:44px;display:flex;align-items:center;gap:12px;padding:6px 10px;background:#17365d;color:#fff;font:700 14px Arial,sans-serif}
.demo-status{margin-right:auto;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.demo-new-tab,.demo-exit{border:1px solid #fff;border-radius:4px;padding:5px 9px;background:#fff;color:#17365d;font:700 13px Arial,sans-serif;text-decoration:none;white-space:nowrap}
.demo-exit{cursor:pointer}
.demo iframe{display:block;box-sizing:border-box;width:100%;height:calc(100% - 44px);border:0;background:#fff}
@media (max-width:900px){.demo-toolbar{height:38px}.demo iframe{height:calc(100% - 38px)}.demo-status{display:none}}
"""


RUNTIME_JS = """(() => {
  function closeDemo(demo) {
    if (!demo) return;
    const frame = demo.querySelector(".demo-frame");
    const activate = demo.querySelector(".demo-activate");
    frame.hidden = true;
    activate.hidden = false;
    activate.focus({ preventScroll: true });
    if (window.Reveal) Reveal.focus();
  }

  document.querySelectorAll(".demo").forEach((demo) => {
    const activate = demo.querySelector(".demo-activate");
    const frame = demo.querySelector(".demo-frame");
    const iframe = demo.querySelector("iframe");
    const status = demo.querySelector(".demo-status");
    activate.addEventListener("click", () => {
      activate.hidden = true;
      frame.hidden = false;
      status.textContent = "Loading interactive demo...";
      if (!iframe.hasAttribute("src")) iframe.src = demo.dataset.demoUrl;
      iframe.focus();
    });
    iframe.addEventListener("load", () => {
      status.textContent = "Interactive demo ready";
    });
    demo.querySelector(".demo-exit").addEventListener("click", () => closeDemo(demo));
  });

  if (window.Reveal) {
    Reveal.on("slidechanged", (event) => {
      if (event.previousSlide) closeDemo(event.previousSlide.querySelector(".demo"));
    });
  }
})();
"""


def write_runtime_assets(output: Path, local_demo: Path) -> None:
    assets = output / "assets"
    demos = output / "demos"
    assets.mkdir(parents=True, exist_ok=True)
    demos.mkdir(parents=True, exist_ok=True)
    (assets / "presentation.css").write_text(RUNTIME_CSS, encoding="ascii")
    (assets / "presentation.js").write_text(RUNTIME_JS, encoding="ascii")
    shutil.copyfile(
        local_demo,
        demos / "IRF4_Fibroblast_top100_log2fc0p5_network.html",
    )


def _png_dimensions(path: Path) -> tuple[int, int]:
    with path.open("rb") as stream:
        header = stream.read(24)
    if len(header) < 24 or not header.startswith(bytes.fromhex("89504e470d0a1a0a")):
        raise ValueError(f"Rendered slide is not a PNG: {path}")
    return struct.unpack(">II", header[16:24])


def validate_rendered_slides(slides: Path, expected_count: int = 59) -> None:
    rendered = sorted(slides.glob("slide-*.png"))
    if len(rendered) != expected_count:
        raise ValueError(f"Expected {expected_count} rendered slides, found {len(rendered)}")
    expected_names = [f"slide-{number:03d}.png" for number in range(1, expected_count + 1)]
    names = [path.name for path in rendered]
    if names != expected_names:
        raise ValueError("Rendered slide filenames are not contiguous and ordered")
    for path in rendered:
        dimensions = _png_dimensions(path)
        if dimensions != (3840, 2160):
            raise ValueError(f"Rendered slide must be 3840x2160: {path} is {dimensions[0]}x{dimensions[1]}")


def copy_reveal_vendor(reveal_root: Path, output: Path) -> None:
    files = {
        reveal_root / "dist" / "reveal.css": output / "vendor" / "reveal" / "reveal.css",
        reveal_root / "dist" / "reveal.js": output / "vendor" / "reveal" / "reveal.js",
        reveal_root / "dist" / "plugin" / "notes.js": output
        / "vendor"
        / "reveal"
        / "plugin"
        / "notes"
        / "notes.js",
        reveal_root / "LICENSE": output / "vendor" / "reveal" / "LICENSE",
    }
    for source, destination in files.items():
        if not source.is_file():
            raise FileNotFoundError(f"Required Reveal.js runtime file is missing: {source}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
    notes_destination = output / "vendor" / "reveal" / "plugin" / "notes" / "notes.js"
    notes_text = notes_destination.read_text(encoding="utf-8")
    notes_destination.write_text(
        notes_text.encode("ascii", "xmlcharrefreplace").decode("ascii"),
        encoding="ascii",
    )


def render_slides(
    source: Path,
    destination: Path,
    expected_count: int = 59,
    soffice: str = "soffice",
    pdftoppm: str = "pdftoppm",
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="craftgrn-slides-") as temporary_directory:
        temporary = Path(temporary_directory)
        profile = temporary / "libreoffice-profile"
        profile.mkdir()
        subprocess.run(
            [
                soffice,
                "--headless",
                "--nologo",
                "--nodefault",
                "--nofirststartwizard",
                f"-env:UserInstallation={profile.resolve().as_uri()}",
                "--convert-to",
                "pdf",
                "--outdir",
                str(temporary),
                str(source.resolve()),
            ],
            check=True,
        )
        pdf = temporary / f"{source.stem}.pdf"
        if not pdf.is_file():
            raise RuntimeError(f"LibreOffice did not create the expected PDF: {pdf}")
        prefix = temporary / "slide"
        subprocess.run(
            [
                pdftoppm,
                "-png",
                "-scale-to-x",
                "3840",
                "-scale-to-y",
                "2160",
                str(pdf),
                str(prefix),
            ],
            check=True,
        )
        _publish_rendered_slides(temporary, destination, expected_count)
    validate_rendered_slides(destination, expected_count=expected_count)


def _host_renderer_provenance(
    soffice: str = "soffice",
    pdftoppm: str = "pdftoppm",
) -> dict[str, str]:
    libreoffice = subprocess.run(
        [soffice, "--version"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    poppler = subprocess.run(
        [pdftoppm, "-v"],
        check=True,
        capture_output=True,
        text=True,
    ).stderr.strip().splitlines()[0]
    return {"mode": "host", "libreoffice": libreoffice, "poppler_utils": poppler}


def _reused_renderer_provenance(output: Path) -> dict[str, str]:
    manifest_path = output / "manifest.json"
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="ascii"))
        provenance = manifest.get("renderer")
        slides = manifest.get("slides")
        hashes_match = isinstance(slides, list) and all(
            isinstance(slide, dict)
            and isinstance(slide.get("number"), int)
            and isinstance(slide.get("image_sha256"), str)
            and hashlib.sha256(
                (
                    output
                    / "assets"
                    / "slides"
                    / f"slide-{slide['number']:03d}.png"
                ).read_bytes()
            ).hexdigest()
            == slide["image_sha256"]
            for slide in slides
        )
        if isinstance(provenance, dict) and provenance and hashes_match:
            return {**provenance, "reused": "true"}
    return {"mode": "reused", "status": "unknown"}


def _publish_rendered_slides(temporary: Path, destination: Path, expected_count: int) -> None:
    rendered = sorted(
        temporary.glob("slide-*.png"),
        key=lambda path: int(path.stem.rsplit("-", 1)[-1]),
    )
    if len(rendered) != expected_count:
        raise ValueError(f"Expected {expected_count} rendered slides, found {len(rendered)}")
    for path in rendered:
        width, height = _png_dimensions(path)
        if (width, height) != (3840, 2160):
            raise ValueError(f"Rendered slide must be 3840x2160: {path} is {width}x{height}")
    destination.mkdir(parents=True, exist_ok=True)
    for old_path in destination.glob("slide-*.png"):
        old_path.unlink()
    for number, path in enumerate(rendered, start=1):
        shutil.copyfile(path, destination / f"slide-{number:03d}.png")


def render_slides_docker(
    source: Path,
    destination: Path,
    expected_count: int = 59,
    docker: str = "docker",
    image: str = RENDERER_IMAGE,
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="craftgrn-slides-") as temporary_directory:
        temporary = Path(temporary_directory)
        created = subprocess.run(
            [docker, "create", image, source.name],
            check=True,
            capture_output=True,
            text=True,
        )
        container = created.stdout.strip()
        if not container:
            raise RuntimeError("Docker did not return a renderer container ID")
        try:
            subprocess.run(
                [docker, "cp", str(source.resolve()), f"{container}:/source/{source.name}"],
                check=True,
            )
            subprocess.run([docker, "start", "-a", container], check=True)
            subprocess.run([docker, "cp", f"{container}:/output/.", str(temporary)], check=True)
        finally:
            subprocess.run(
                [docker, "rm", "-f", container],
                check=False,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        _publish_rendered_slides(temporary, destination, expected_count)
    validate_rendered_slides(destination, expected_count=expected_count)


def build_deck(
    source: Path,
    output: Path,
    local_demo: Path,
    reveal_root: Path,
    expected_slide_count: int = 59,
    expected_source_sha256: str | None = None,
    render: bool = True,
    use_docker: bool = True,
) -> None:
    source_sha256 = hashlib.sha256(source.read_bytes()).hexdigest()
    if expected_source_sha256 is not None and source_sha256 != expected_source_sha256:
        raise ValueError(
            "PowerPoint source checksum does not match the approved deck: "
            f"expected {expected_source_sha256}, found {source_sha256}"
        )
    presentation = parse_presentation(source)
    if len(presentation.slides) != expected_slide_count:
        raise ValueError(f"Expected {expected_slide_count} slides, found {len(presentation.slides)}")
    slides = output / "assets" / "slides"
    if render:
        if use_docker:
            render_slides_docker(source, slides, expected_count=expected_slide_count)
            renderer_provenance = RENDERER_PROVENANCE
        else:
            render_slides(source, slides, expected_count=expected_slide_count)
            renderer_provenance = _host_renderer_provenance()
    else:
        validate_rendered_slides(slides, expected_count=expected_slide_count)
        renderer_provenance = _reused_renderer_provenance(output)
    generate_site(
        presentation,
        source,
        output,
        expected_slide_count=expected_slide_count,
        renderer_provenance=renderer_provenance,
    )
    write_runtime_assets(output, local_demo)
    copy_reveal_vendor(reveal_root, output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    presentation_root = Path(__file__).resolve().parents[1]
    parser.add_argument(
        "--source",
        type=Path,
        default=presentation_root / "20260813_thesis_presentation.pptx",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=presentation_root / "reveal",
    )
    parser.add_argument(
        "--local-demo",
        type=Path,
        default=presentation_root / "IRF4_Fibroblast_top100_log2fc0p5_network.html",
    )
    parser.add_argument(
        "--reveal-root",
        type=Path,
        default=presentation_root / "node_modules" / "reveal.js",
    )
    parser.add_argument("--skip-render", action="store_true")
    parser.add_argument("--host-renderer", action="store_true")
    args = parser.parse_args()
    build_deck(
        args.source,
        args.output,
        args.local_demo,
        args.reveal_root,
        expected_source_sha256=EXPECTED_SOURCE_SHA256,
        render=not args.skip_render,
        use_docker=not args.host_renderer,
    )


if __name__ == "__main__":
    main()
