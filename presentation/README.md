# Thesis presentation

This directory contains the generated Reveal.js version of the 2026 thesis
presentation. The PowerPoint source remains local and is ignored by Git.

## Build

Install Node dependencies and build the pinned rendering container:

```sh
npm ci
npm run renderer:build
npm run build
```

The build validates the approved PowerPoint SHA-256 checksum, extracts all
speaker notes, renders 68 slides at 3840x2160, and writes the public site to
`reveal/`. Use `npm run build:site` to regenerate HTML and runtime assets from
an existing validated slide render.

The renderer pins its Ubuntu base image by digest, uses Canonical's immutable
Ubuntu archive snapshot, and pins the LibreOffice, Poppler, and font package
versions used for rasterization. The generated
manifest records renderer provenance and a SHA-256 checksum for every slide
image. The renderer substitutes Liberation Sans for Arial, Carlito for Calibri,
Caladea for Cambria, and Liberation Sans for Google Sans, Segoe UI, and
Tahoma. These metric-compatible substitutions make headless rendering
reproducible on systems without the original proprietary fonts. To use a
locally installed LibreOffice and Poppler instead, pass `--host-renderer`; the
manifest records the host tool versions rather than Docker provenance.

Serve the repository root to test the same `/presentation/reveal/` path used
by the browser tests:

```sh
python3 -m http.server 8000 --directory ..
npm run test:e2e
```

In the presentation, press `S` to open Reveal.js speaker view. Slides 14 and 22
expose click-to-activate interactive demonstrations; every other slide is a
static 3840x2160 image.
