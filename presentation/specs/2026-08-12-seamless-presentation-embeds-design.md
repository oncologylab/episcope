# Seamless Presentation Embeds Design

## Scope

Refine the existing Reveal.js embeds on slides 14 and 43 without changing the
PowerPoint-derived slide images, slide numbering, speaker notes, or slide 20.

## Slide 14

The static PowerPoint screenshot remains visible before activation. Activating
the demonstration replaces exactly the screenshot region with the complete
remote differential-footprint dashboard, including its configuration, export,
selection, and plotting controls.

The embedded page uses a larger virtual viewport and is scaled down to match the
layout and density of the screenshot already present on the slide. The Reveal
embed must not add a persistent toolbar, border, shadow, title, or background.
Small translucent `Open` and `Exit` controls float over the upper-right corner
and become prominent on pointer hover or keyboard focus. The existing new-tab
fallback remains available.

The remote URL remains the source of the interactive dashboard. The deck does
not copy or alter the remote report and retains the existing sandbox boundary.

## Slide 43

The local network report gains an explicit canvas-only embed mode. In this mode
the page retains its network data and interaction code but hides the report
heading, metadata, controls, statistics, outer padding, borders, and shadow.
The `#canvas` element fills the iframe viewport edge-to-edge with a transparent
or white background matching the slide.

The slide uses this mode and places the canvas in the existing screenshot
region. Network dragging, node selection, tooltips, and SVG rendering continue
to work. The same unobtrusive floating `Open` and `Exit` controls used on slide
14 provide recovery without making the embed look like a separate webpage.

## Shared Behavior

- Keep click-to-activate loading so remote content is not requested early.
- Keep direct new-tab fallbacks and keyboard-focus visibility.
- Close an active embed when leaving its slide and return focus to Reveal.
- Preserve iframe sandboxing and do not grant same-origin access to remote code.
- Preserve all 59 speaker-note blocks and the Reveal speaker console.
- Leave slide 20 behavior unchanged.

## Implementation Boundaries

The presentation builder remains the source of generated HTML, CSS, and
JavaScript. Per-demo configuration describes the display mode and slide 14
scale. The local network source file implements canvas-only mode through a URL
query parameter or fragment-derived class; the builder copies that source into
the generated site as it does today.

Generated presentation assets are rebuilt from these sources and committed so
GitHub Pages can publish them without running the presentation renderer.

## Validation

Automated tests must verify:

- Slide 14 replaces the screenshot region without the persistent toolbar and
  uses the configured scaled viewport while the complete mocked dashboard stays
  interactive.
- Slide 43 hides all report chrome, makes the canvas fill the iframe viewport,
  and retains graph interaction.
- Floating recovery controls are keyboard accessible.
- Slide 20, fallback behavior, iframe isolation, slide navigation, and speaker
  notes continue to pass their existing tests.

Final validation includes desktop Playwright screenshots for slides 14 and 43,
the full unit and browser suites, source/package checks proportional to the
changed files, GitHub Pages deployment, and a browser check of the live deck.
