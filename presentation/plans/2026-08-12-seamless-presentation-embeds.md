# Seamless Presentation Embeds Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make slide 14's complete live dashboard match its PowerPoint screenshot and make slide 43 display only its interactive network canvas.

**Architecture:** Keep the presentation builder as the source of generated markup and runtime assets. Add per-demo presentation modes so slide 14 receives a 0.5-scale virtual viewport and slides 14/43 use overlay controls instead of a toolbar; add a query-driven canvas-only mode to the bundled slide-43 report.

**Tech Stack:** Python 3 builder and unittest, Reveal.js 6, HTML/CSS/JavaScript, Playwright 1.62.1, GitHub Pages.

## Global Constraints

- Preserve all 59 PowerPoint-derived 3840x2160 slide images and speaker notes.
- Retain the slide 14 remote URL and iframe sandbox without `allow-same-origin`.
- Keep click-to-activate loading, direct-link fallback, and keyboard recovery.
- Leave slide 20 behavior unchanged.
- Generate source and runtime text as ASCII only.
- Do not stage or modify unrelated working-tree changes.

---

### Task 1: Per-demo seamless framing

**Files:**
- Modify: `presentation/tests/test_build_presentation.py`
- Modify: `presentation/tests/presentation.spec.js`
- Modify: `presentation/tools/build_presentation.py`
- Regenerate: `presentation/reveal/index.html`
- Regenerate: `presentation/reveal/assets/presentation.css`

**Interfaces:**
- Consumes: `DEMOS`, `_demo_markup(slide_number)`, and `RUNTIME_CSS`.
- Produces: `data-demo-mode`, `--demo-scale`, `.demo-seamless`, and overlay toolbar behavior used by slides 14 and 43.

- [ ] **Step 1: Add failing unit and browser assertions**

Assert that slide 14 markup includes `data-demo-mode="scaled-dashboard"`,
`--demo-scale:0.5`, and a seamless class; assert slide 20 retains its standard
mode. In Playwright, activate slide 14 and assert its toolbar is an absolute
overlay, its frame has no border, and its iframe rendered width is twice the
visible frame width.

- [ ] **Step 2: Run the focused tests and verify the new assertions fail**

Run:

```bash
cd presentation
python3 -m unittest tests.test_build_presentation.PresentationBuilderTests.test_generates_site_with_notes_manifest_and_demo_contract -v
npx playwright test tests/presentation.spec.js --grep "remote reports"
```

Expected: failures because seamless demo modes and scaling are absent.

- [ ] **Step 3: Implement per-demo modes and seamless CSS**

Add `mode` and `scale` keys to the slide 14 and 43 demo definitions. Generate
mode-specific classes and custom properties in `_demo_markup()`. Keep the
existing toolbar markup for accessibility, but for `.demo-seamless` position it
as a small translucent upper-right overlay, visually hide the status text,
remove the frame border/shadow, clip overflow, and render the iframe as:

```css
width:calc(100% / var(--demo-scale));
height:calc(100% / var(--demo-scale));
transform:scale(var(--demo-scale));
transform-origin:top left;
```

Keep standard toolbar sizing for slide 20.

- [ ] **Step 4: Rebuild generated assets and run the focused tests**

Run:

```bash
cd presentation
npm run build:site
npm run test:unit
npx playwright test tests/presentation.spec.js --grep "remote reports"
```

Expected: all commands exit successfully.

- [ ] **Step 5: Commit the framing change**

Stage only the Task 1 source, tests, and generated presentation assets, then
commit with message `Refine interactive presentation framing`.

### Task 2: Canvas-only slide 43 mode

**Files:**
- Modify: `presentation/tests/presentation.spec.js`
- Modify: `presentation/IRF4_Fibroblast_top100_log2fc0p5_network.html`
- Regenerate: `presentation/reveal/demos/IRF4_Fibroblast_top100_log2fc0p5_network.html`
- Regenerate: `presentation/reveal/index.html`

**Interfaces:**
- Consumes: slide 43 URL `demos/IRF4_Fibroblast_top100_log2fc0p5_network.html?embed=canvas`.
- Produces: the `canvas-embed` document class and full-viewport `#canvas` layout.

- [ ] **Step 1: Add a failing browser assertion for canvas-only layout**

After activating slide 43, assert the iframe URL contains `embed=canvas`, the
`.top` and `#stats` elements are hidden, and the `#canvas` bounding box matches
the iframe viewport within one CSS pixel. Keep the existing graph SVG, layout
selection, and navigation assertions.

- [ ] **Step 2: Run the focused test and verify it fails for missing canvas mode**

Run:

```bash
cd presentation
npx playwright test tests/presentation.spec.js --grep "local IRF4 network"
```

Expected: failure because the current page shows its title, controls, and stats.

- [ ] **Step 3: Implement query-driven canvas-only layout**

Change the slide 43 demo URL to append `?embed=canvas`. In the local report,
detect `new URLSearchParams(location.search).get("embed") === "canvas"` and add
`canvas-embed` to the root element. Add CSS that makes `html`, `body`, `.wrap`,
and `#canvas` fill the viewport, hides `.top` and `#stats`, removes padding,
border, shadow, aspect-ratio, and maximum-height constraints, and keeps the SVG
and tooltip interaction unchanged.

- [ ] **Step 4: Rebuild and run the focused test**

Run:

```bash
cd presentation
npm run build:site
npx playwright test tests/presentation.spec.js --grep "local IRF4 network"
```

Expected: the focused Playwright test passes.

- [ ] **Step 5: Commit the canvas-only change**

Stage only the Task 2 source, tests, and generated presentation assets, then
commit with message `Show canvas-only network on slide 43`.

### Task 3: Full verification and publication

**Files:**
- Verify: `presentation/reveal/`
- Verify: `.github/workflows/pkgdown.yaml`
- Verify: `tests/testthat/test-config-pkgdown.R`

**Interfaces:**
- Consumes: the complete generated presentation tree.
- Produces: a verified GitHub Pages presentation at `/craftgrn/presentation/`.

- [ ] **Step 1: Run all presentation tests**

Run:

```bash
cd presentation
npm run test:unit
npm run test:e2e
```

Expected: 9 unit tests and all Playwright tests pass.

- [ ] **Step 2: Capture and inspect desktop screenshots**

Use Playwright at 1920x1080 to activate slides 14 and 43, wait for their live
content, and save screenshots under `/tmp`. Confirm slide 14 aligns with its
static screenshot density and slide 43 contains only the network canvas plus
the unobtrusive recovery controls.

- [ ] **Step 3: Run repository checks proportional to the change**

Run:

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-config-pkgdown.R")'
git diff --check
find presentation/tools presentation/tests presentation/reveal -type f -not -path '*/vendor/*' -print0 | xargs -0 grep -Il $'[^\x00-\x7F]' || true
```

Expected: the pkgdown configuration test passes, no whitespace errors are
reported, and the ASCII audit lists no changed source/runtime file.

- [ ] **Step 4: Review, commit any final generated changes, and push main**

Confirm `HEAD` is on `main`, inspect the complete presentation diff, stage only
presentation-owned files, commit any remaining generated changes, and push
`main` to `origin`.

- [ ] **Step 5: Monitor GitHub Pages and validate the live deck**

Wait for the pushed pkgdown workflow and deployment jobs to succeed. Then run a
live Playwright check that verifies 59 slides, 59 note blocks, slide 14's scaled
dashboard, slide 43's canvas-only layout and graph interaction, and speaker
view notes.
