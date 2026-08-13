import importlib.util
import json
import struct
import os
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "presentation" / "tools" / "build_presentation.py"


def load_builder():
    spec = importlib.util.spec_from_file_location("build_presentation", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def write_fixture(path, height=6858000):
    presentation = """<?xml version="1.0" encoding="UTF-8"?>
<p:presentation xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main"
 xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
 <p:sldIdLst><p:sldId id="256" r:id="rId2"/><p:sldId id="257" r:id="rId1"/></p:sldIdLst>
 <p:sldSz cx="12192000" cy="{height}"/>
</p:presentation>""".format(height=height)
    presentation_rels = """<?xml version="1.0" encoding="UTF-8"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
 <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slide" Target="slides/slide1.xml"/>
 <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slide" Target="slides/slide2.xml"/>
</Relationships>"""
    slide = """<?xml version="1.0" encoding="UTF-8"?>
<p:sld xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main"
 xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">
 <p:cSld><p:spTree><p:sp><p:nvSpPr><p:nvPr><p:ph type="title"/></p:nvPr></p:nvSpPr>
 <p:txBody><a:p><a:r><a:t>{title}</a:t></a:r></a:p></p:txBody></p:sp></p:spTree></p:cSld>
</p:sld>"""
    slide_rels = """<?xml version="1.0" encoding="UTF-8"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
 <Relationship Id="rId9" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/notesSlide" Target="../notesSlides/notesSlide2.xml"/>
</Relationships>"""
    notes = """<?xml version="1.0" encoding="UTF-8"?>
<p:notes xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main"
 xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main">
 <p:cSld><p:spTree>
  <p:sp><p:nvSpPr><p:nvPr><p:ph type="body"/></p:nvPr></p:nvSpPr><p:txBody>
   <a:p><a:r><a:t>First paragraph.</a:t></a:r></a:p>
   <a:p><a:r><a:rPr b="1"/><a:t>Bold</a:t></a:r><a:r><a:t> and </a:t></a:r><a:r><a:rPr i="1"/><a:t>italic</a:t></a:r><a:r><a:t> &#x2014; exact.</a:t></a:r></a:p>
  </p:txBody></p:sp>
  <p:sp><p:nvSpPr><p:nvPr><p:ph type="sldNum"/></p:nvPr></p:nvSpPr><p:txBody><a:p><a:r><a:t>2</a:t></a:r></a:p></p:txBody></p:sp>
 </p:spTree></p:cSld>
</p:notes>"""
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("ppt/presentation.xml", presentation)
        archive.writestr("ppt/_rels/presentation.xml.rels", presentation_rels)
        archive.writestr("ppt/slides/slide1.xml", slide.format(title="Physical slide one"))
        archive.writestr("ppt/slides/slide2.xml", slide.format(title="First in show"))
        archive.writestr("ppt/slides/_rels/slide2.xml.rels", slide_rels)
        archive.writestr("ppt/notesSlides/notesSlide2.xml", notes)


def write_png_header(path, width, height):
    path.write_bytes(
        bytes.fromhex("89504e470d0a1a0a")
        + b"\x00\x00\x00\x0dIHDR"
        + struct.pack(">II", width, height)
        + b"\x08\x06\x00\x00\x00"
    )


class ParsePresentationTests(unittest.TestCase):
    def test_uses_presentation_relationship_order_and_preserves_note_formatting(self):
        builder = load_builder()
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "fixture.pptx"
            write_fixture(source)
            parsed = builder.parse_presentation(source)

        self.assertEqual((parsed.width, parsed.height), (12192000, 6858000))
        self.assertEqual([slide.source_number for slide in parsed.slides], [2, 1])
        self.assertEqual([slide.title for slide in parsed.slides], ["First in show", "Physical slide one"])
        self.assertEqual(
            parsed.slides[0].notes_text,
            "First paragraph.\nBold and italic \u2014 exact.",
        )
        self.assertEqual(
            parsed.slides[0].notes_html,
            "<p>First paragraph.</p><p><strong>Bold</strong> and <em>italic</em> &#8212; exact.</p>",
        )
        self.assertEqual(parsed.slides[1].notes_text, "")
        self.assertEqual(parsed.slides[1].notes_html, "")

    def test_rejects_non_widescreen_source(self):
        builder = load_builder()
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "fixture.pptx"
            write_fixture(source, height=9144000)
            with self.assertRaisesRegex(ValueError, "16:9"):
                builder.parse_presentation(source)


class SiteGenerationTests(unittest.TestCase):
    def test_generates_reveal_sections_manifest_and_exact_demo_mappings(self):
        builder = load_builder()
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / "fixture.pptx"
            output = directory / "site"
            write_fixture(source)
            parsed = builder.parse_presentation(source)
            slides = output / "assets" / "slides"
            slides.mkdir(parents=True)
            write_png_header(slides / "slide-001.png", 3840, 2160)
            write_png_header(slides / "slide-002.png", 3840, 2160)
            builder.generate_site(
                parsed,
                source,
                output,
                expected_slide_count=2,
                renderer_provenance=builder.RENDERER_PROVENANCE,
            )

            html = (output / "index.html").read_text(encoding="ascii")
            manifest = json.loads((output / "manifest.json").read_text(encoding="ascii"))

        self.assertEqual(html.count("<section "), 2)
        self.assertIn('<aside class="notes"><p>First paragraph.</p>', html)
        self.assertIn('slideNumber: "c/t"', html)
        self.assertEqual(manifest["slide_count"], 2)
        self.assertEqual(manifest["slide_size_emu"], [12192000, 6858000])
        self.assertRegex(manifest["source_sha256"], r"^[0-9a-f]{64}$")
        self.assertEqual(manifest["renderer"], builder.RENDERER_PROVENANCE)
        self.assertRegex(manifest["slides"][0]["image_sha256"], r"^[0-9a-f]{64}$")
        self.assertIn(
            'sandbox="allow-scripts allow-forms allow-downloads allow-popups allow-popups-to-escape-sandbox"',
            builder._demo_markup(14),
        )
        self.assertNotIn("allow-same-origin", builder._demo_markup(14))
        self.assertIn('class="demo demo-seamless demo-scaled-dashboard"', builder._demo_markup(14))
        self.assertIn('data-demo-mode="scaled-dashboard"', builder._demo_markup(14))
        self.assertIn("--demo-scale:0.5", builder._demo_markup(14))
        self.assertIn('data-demo-mode="standard"', builder._demo_markup(20))
        self.assertNotIn("demo-seamless", builder._demo_markup(20))
        self.assertEqual(
            builder.DEMOS[14]["url"],
            "https://oncologylab.github.io/fp-tools/demos/reports/diff_footprints_K562_HepG2.html",
        )
        self.assertEqual(
            builder.DEMOS[20]["url"],
            "https://oncologylab.github.io/fp-tools/demos/gui/fp-tools-gui-static-demo.html#home",
        )
        self.assertEqual(
            builder.DEMOS[43]["url"],
            "demos/IRF4_Fibroblast_top100_log2fc0p5_network.html?embed=canvas",
        )

    def test_writes_runtime_assets_and_copies_local_demo_byte_for_byte(self):
        builder = load_builder()
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            output = directory / "site"
            local_demo = directory / "network.html"
            local_demo.write_bytes(b"<!doctype html><title>network</title>\n")
            builder.write_runtime_assets(output, local_demo)

            javascript = (output / "assets" / "presentation.js").read_text(encoding="ascii")
            stylesheet = (output / "assets" / "presentation.css").read_text(encoding="ascii")
            copied_demo = output / "demos" / "IRF4_Fibroblast_top100_log2fc0p5_network.html"
            copied_demo_bytes = copied_demo.read_bytes()

        self.assertIn("iframe.src = demo.dataset.demoUrl", javascript)
        self.assertIn("Reveal.on(\"slidechanged\"", javascript)
        self.assertIn("--demo-left", stylesheet)
        self.assertEqual(copied_demo_bytes, b"<!doctype html><title>network</title>\n")

    def test_render_validation_requires_every_4k_slide(self):
        builder = load_builder()
        with tempfile.TemporaryDirectory() as directory:
            slides = Path(directory)
            write_png_header(slides / "slide-001.png", 3840, 2160)
            write_png_header(slides / "slide-002.png", 3840, 2160)
            builder.validate_rendered_slides(slides, expected_count=2)
            write_png_header(slides / "slide-002.png", 1920, 1080)
            with self.assertRaisesRegex(ValueError, "3840x2160"):
                builder.validate_rendered_slides(slides, expected_count=2)
            (slides / "slide-002.png").unlink()
            with self.assertRaisesRegex(ValueError, "Expected 2 rendered slides"):
                builder.validate_rendered_slides(slides, expected_count=2)

    def test_copies_only_required_reveal_runtime_files(self):
        builder = load_builder()
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            reveal = directory / "node_modules" / "reveal.js"
            (reveal / "dist" / "plugin").mkdir(parents=True)
            (reveal / "dist" / "reveal.css").write_text("reveal css\n", encoding="ascii")
            (reveal / "dist" / "reveal.js").write_text("reveal js\n", encoding="ascii")
            (reveal / "dist" / "plugin" / "notes.js").write_text("notes \u2013 js\n", encoding="utf-8")
            (reveal / "LICENSE").write_text("MIT license\n", encoding="ascii")
            (reveal / "dist" / "unneeded.js").write_text("unused\n", encoding="ascii")
            output = directory / "site"
            builder.copy_reveal_vendor(reveal, output)

            vendor = output / "vendor" / "reveal"
            self.assertEqual((vendor / "reveal.css").read_text(), "reveal css\n")
            self.assertEqual((vendor / "reveal.js").read_text(), "reveal js\n")
            notes_bytes = (vendor / "plugin" / "notes" / "notes.js").read_bytes()
            self.assertEqual(notes_bytes, b"notes &#8211; js\n")
            self.assertEqual((vendor / "LICENSE").read_text(), "MIT license\n")
            self.assertFalse((vendor / "unneeded.js").exists())

    def test_renderer_normalizes_poppler_pages_to_contiguous_4k_names(self):
        builder = load_builder()
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / "deck.pptx"
            source.write_bytes(b"pptx")
            bin_dir = directory / "bin"
            bin_dir.mkdir()
            fake_tool = bin_dir / "fake-render-tool"
            fake_tool.write_text(
                """#!/usr/bin/env python3
import pathlib, struct, sys
name = pathlib.Path(sys.argv[0]).name
if name == 'soffice':
    out = pathlib.Path(sys.argv[sys.argv.index('--outdir') + 1])
    out.joinpath('deck.pdf').write_bytes(b'%PDF-fake')
else:
    prefix = pathlib.Path(sys.argv[-1])
    header = bytes.fromhex('89504e470d0a1a0a') + b'\\0\\0\\0\\rIHDR' + struct.pack('>II', 3840, 2160) + b'\\x08\\x06\\0\\0\\0'
    prefix.with_name(prefix.name + '-1.png').write_bytes(header)
    prefix.with_name(prefix.name + '-2.png').write_bytes(header)
""",
                encoding="ascii",
            )
            fake_tool.chmod(0o755)
            os.symlink(fake_tool, bin_dir / "soffice")
            os.symlink(fake_tool, bin_dir / "pdftoppm")
            slides = directory / "slides"
            builder.render_slides(
                source,
                slides,
                expected_count=2,
                soffice=str(bin_dir / "soffice"),
                pdftoppm=str(bin_dir / "pdftoppm"),
            )

            self.assertEqual(
                sorted(path.name for path in slides.glob("*.png")),
                ["slide-001.png", "slide-002.png"],
            )
            builder.validate_rendered_slides(slides, expected_count=2)

    def test_build_deck_assembles_existing_renders_and_rejects_wrong_source(self):
        builder = load_builder()
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / "fixture.pptx"
            write_fixture(source)
            local_demo = directory / "network.html"
            local_demo.write_text("<!doctype html>network\n", encoding="ascii")
            reveal = directory / "reveal.js"
            (reveal / "dist" / "plugin").mkdir(parents=True)
            for path in (
                reveal / "dist" / "reveal.css",
                reveal / "dist" / "reveal.js",
                reveal / "dist" / "plugin" / "notes.js",
                reveal / "LICENSE",
            ):
                path.write_text(path.name + "\n", encoding="ascii")
            output = directory / "site"
            slides = output / "assets" / "slides"
            slides.mkdir(parents=True)
            write_png_header(slides / "slide-001.png", 3840, 2160)
            write_png_header(slides / "slide-002.png", 3840, 2160)

            builder.build_deck(
                source,
                output,
                local_demo,
                reveal,
                expected_slide_count=2,
                render=False,
            )
            self.assertTrue((output / "index.html").is_file())
            self.assertTrue((output / "assets" / "presentation.js").is_file())
            self.assertTrue((output / "vendor" / "reveal" / "reveal.js").is_file())
            manifest = json.loads((output / "manifest.json").read_text(encoding="ascii"))
            self.assertEqual(manifest["renderer"], {"mode": "reused", "status": "unknown"})

            manifest["renderer"] = builder.RENDERER_PROVENANCE
            (output / "manifest.json").write_text(json.dumps(manifest), encoding="ascii")
            write_png_header(slides / "slide-001.png", 3840, 2160)
            with (slides / "slide-001.png").open("ab") as stream:
                stream.write(b"changed")
            self.assertEqual(
                builder._reused_renderer_provenance(output),
                {"mode": "reused", "status": "unknown"},
            )

            with self.assertRaisesRegex(ValueError, "checksum"):
                builder.build_deck(
                    source,
                    output,
                    local_demo,
                    reveal,
                    expected_slide_count=2,
                    expected_source_sha256="0" * 64,
                    render=False,
                )

    def test_docker_renderer_copies_source_through_daemon_and_produces_4k_slides(self):
        builder = load_builder()
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / "deck.pptx"
            source.write_bytes(b"pptx")
            fake_docker = directory / "docker"
            fake_docker.write_text(
                """#!/usr/bin/env python3
import pathlib, shutil, struct, sys
root = pathlib.Path(__file__).with_name('docker-state')
root.mkdir(exist_ok=True)
log = pathlib.Path(__file__).with_name('docker.log')
with log.open('a') as stream: stream.write(' '.join(sys.argv[1:]) + '\\n')
command = sys.argv[1]
if command == 'create':
    (root / 'source').mkdir(exist_ok=True)
    (root / 'output').mkdir(exist_ok=True)
    print('fake-container')
elif command == 'cp':
    source, destination = sys.argv[2:4]
    if ':' in destination:
        shutil.copyfile(source, root / 'source' / pathlib.Path(destination).name)
    else:
        pathlib.Path(destination).mkdir(parents=True, exist_ok=True)
        for path in (root / 'output').iterdir(): shutil.copyfile(path, pathlib.Path(destination) / path.name)
elif command == 'start':
    assert (root / 'source' / 'deck.pptx').read_bytes() == b'pptx'
    header = bytes.fromhex('89504e470d0a1a0a') + b'\\0\\0\\0\\rIHDR' + struct.pack('>II', 3840, 2160) + b'\\x08\\x06\\0\\0\\0'
    (root / 'output' / 'slide-1.png').write_bytes(header)
    (root / 'output' / 'slide-2.png').write_bytes(header)
elif command == 'rm':
    pass
else:
    raise SystemExit('unexpected docker command: ' + command)
""",
                encoding="ascii",
            )
            fake_docker.chmod(0o755)
            slides = directory / "slides"
            builder.render_slides_docker(
                source,
                slides,
                expected_count=2,
                docker=str(fake_docker),
                image="craftgrn-test-renderer",
            )
            self.assertEqual(
                sorted(path.name for path in slides.glob("*.png")),
                ["slide-001.png", "slide-002.png"],
            )
            docker_log = (directory / "docker.log").read_text(encoding="ascii")
            self.assertIn("create", docker_log)
            self.assertIn("cp", docker_log)
            self.assertIn("start -a", docker_log)
            self.assertNotIn(" -v ", f" {docker_log} ")


if __name__ == "__main__":
    unittest.main()
