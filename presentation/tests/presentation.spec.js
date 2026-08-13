const { test, expect } = require("@playwright/test");

async function openDeck(page) {
  await page.goto("./");
  await page.waitForFunction(() => window.Reveal && Reveal.isReady());
}

async function goToSlide(page, number) {
  await page.evaluate((target) => Reveal.slide(target - 1), number);
  await expect(page.locator(`#slide-${number}`)).toHaveClass(/present/);
}

test("loads all 59 slides with matching deep links and numbering", async ({ page }) => {
  await openDeck(page);
  await expect(page.locator(".reveal .slides > section")).toHaveCount(59);
  await expect(page.locator(".slide-number")).toContainText("1 / 59");
  await goToSlide(page, 59);
  await expect(page).toHaveURL(/#\/slide-59$/);
  await expect(page.locator(".slide-number")).toContainText("59 / 59");
});

test("keeps remote reports lazy and provides exit and fallback controls", async ({ page }) => {
  const reportUrl = "https://oncologylab.github.io/fp-tools/demos/reports/diff_footprints_K562_HepG2.html";
  let reportRequests = 0;
  await page.route(reportUrl, async (route) => {
    reportRequests += 1;
    await route.fulfill({ contentType: "text/html", body: "<!doctype html><div id=dashboard>Report ready</div>" });
  });
  await openDeck(page);
  await goToSlide(page, 14);
  expect(reportRequests).toBe(0);
  const slide = page.locator("#slide-14");
  await expect(slide.locator("iframe")).toHaveAttribute(
    "sandbox",
    "allow-scripts allow-forms allow-downloads allow-popups allow-popups-to-escape-sandbox"
  );
  await slide.locator(".demo-activate").click();
  await expect(slide.locator(".demo-frame")).toBeVisible();
  await expect(slide.locator(".demo-fallback")).toBeHidden();
  const seamlessLayout = await slide.locator(".demo-frame").evaluate((frame) => {
    const iframe = frame.querySelector("iframe");
    const toolbar = frame.querySelector(".demo-toolbar");
    const frameStyle = getComputedStyle(frame);
    return {
      frameWidth: frame.clientWidth,
      iframeLayoutWidth: iframe.offsetWidth,
      toolbarPosition: getComputedStyle(toolbar).position,
      frameBorderWidth: frameStyle.borderTopWidth,
      frameBoxShadow: frameStyle.boxShadow,
    };
  });
  expect(seamlessLayout.iframeLayoutWidth).toBeCloseTo(seamlessLayout.frameWidth / 0.9, -1);
  expect(seamlessLayout.toolbarPosition).toBe("absolute");
  expect(seamlessLayout.frameBorderWidth).toBe("0px");
  expect(seamlessLayout.frameBoxShadow).toBe("none");
  await expect(slide.locator(".demo-new-tab")).toHaveAttribute("href", reportUrl);
  await expect(slide.locator("iframe").contentFrame().locator("#dashboard")).toHaveText("Report ready");
  expect(reportRequests).toBe(1);
  await slide.locator(".demo-exit").click();
  await expect(slide.locator(".demo-frame")).toBeHidden();
  await expect(slide.locator(".demo-fallback")).toBeVisible();
  await page.keyboard.press("ArrowRight");
  await expect(page.locator("#slide-15")).toHaveClass(/present/);
});

test("isolates interactive demo scripts from the parent presentation", async ({ page }) => {
  const reportUrl = "https://oncologylab.github.io/fp-tools/demos/reports/diff_footprints_K562_HepG2.html";
  await page.route(reportUrl, async (route) => {
    await route.fulfill({
      contentType: "text/html",
      body: `<!doctype html><output id="probe"></output><script>
        try {
          parent.document.body.dataset.compromised = "yes";
          document.querySelector("#probe").textContent = "parent-accessible";
        } catch (error) {
          document.querySelector("#probe").textContent = "parent-isolated";
        }
      </script>`,
    });
  });
  await openDeck(page);
  await goToSlide(page, 14);
  const slide = page.locator("#slide-14");
  await slide.locator(".demo-activate").click();
  await expect(slide.locator("iframe").contentFrame().locator("#probe")).toHaveText("parent-isolated");
  await expect(page.locator("body")).not.toHaveAttribute("data-compromised", "yes");
});

test("restores the static slide and direct link when a remote demo is blocked", async ({ page }) => {
  const reportUrl = "https://oncologylab.github.io/fp-tools/demos/reports/diff_footprints_K562_HepG2.html";
  await page.route(reportUrl, (route) => route.abort("internetdisconnected"));
  await openDeck(page);
  await goToSlide(page, 14);
  const slide = page.locator("#slide-14");
  await slide.locator(".demo-activate").click();
  await expect(slide.locator(".demo-frame")).toBeVisible();
  await slide.locator(".demo-exit").click();
  await expect(slide.locator(".demo-frame")).toBeHidden();
  await expect(slide.locator(".demo-activate")).toBeVisible();
  await expect(slide.locator(".demo-fallback")).toBeVisible();
  await expect(slide.locator(".demo-fallback")).toHaveAttribute("href", reportUrl);
});

test("loads the GUI and local IRF4 network in their requested slides", async ({ page }) => {
  const guiUrl = "https://oncologylab.github.io/fp-tools/demos/gui/fp-tools-gui-static-demo.html#home";
  await page.route("https://oncologylab.github.io/fp-tools/demos/gui/fp-tools-gui-static-demo.html", async (route) => {
    await route.fulfill({ contentType: "text/html", body: "<!doctype html><main id=main-content>GUI ready</main>" });
  });
  await openDeck(page);
  await goToSlide(page, 20);
  const guiSlide = page.locator("#slide-20");
  await guiSlide.locator(".demo-activate").click();
  await expect(guiSlide.locator("iframe").contentFrame().locator("#main-content")).toHaveText("GUI ready");
  await guiSlide.locator(".demo-exit").click();

  await goToSlide(page, 43);
  const networkSlide = page.locator("#slide-43");
  await networkSlide.locator(".demo-activate").click();
  await expect(networkSlide.locator("iframe")).toHaveAttribute("src", /embed=canvas$/);
  const network = networkSlide.locator("iframe").contentFrame();
  await expect(network.locator(".top")).toBeHidden();
  await expect(network.locator("#stats")).toBeHidden();
  await expect(network.locator("#canvas svg")).toBeVisible();
  const canvasLayout = await network.locator("#canvas").evaluate((canvas) => {
    const bounds = canvas.getBoundingClientRect();
    return {
      left: bounds.left,
      top: bounds.top,
      width: bounds.width,
      height: bounds.height,
      viewportWidth: window.innerWidth,
      viewportHeight: window.innerHeight,
    };
  });
  expect(canvasLayout.left).toBeCloseTo(0, 0);
  expect(canvasLayout.top).toBeCloseTo(0, 0);
  expect(canvasLayout.width).toBeCloseTo(canvasLayout.viewportWidth, 0);
  expect(canvasLayout.height).toBeCloseTo(canvasLayout.viewportHeight, 0);
  await expect(network.locator("#nodeLayer > *").first()).toBeVisible();
  await networkSlide.locator(".demo-exit").click();
  await page.keyboard.press("ArrowRight");
  await expect(page.locator("#slide-44")).toHaveClass(/present/);
});

test("opens speaker view with current notes, next slide, clock, and timer", async ({ page }) => {
  await openDeck(page);
  const popupPromise = page.waitForEvent("popup");
  await page.keyboard.press("s");
  const speaker = await popupPromise;
  await expect(speaker.locator("#current-slide iframe")).toBeVisible();
  await expect(speaker.locator("#upcoming-slide iframe")).toBeVisible();
  await expect(speaker.locator(".speaker-controls-notes .value")).toContainText("Good morning, and thank you for being here.");
  await expect(speaker.locator(".clock-value")).not.toHaveText("0:00 AM");
  await expect(speaker.locator(".timer")).toContainText(":");
});

test("has one aligned notes element per slide, including explicit empty notes", async ({ page }) => {
  await openDeck(page);
  await expect(page.locator("aside.notes")).toHaveCount(59);
  await expect(page.locator("#slide-14 aside.notes")).toContainText("HNF4A and ONECUT2 are stronger in HepG2");
  await expect(page.locator("#slide-20 aside.notes")).toContainText("same reproducible workflow");
  await expect(page.locator("#slide-43 aside.notes")).toContainText("coordinated multi-TF regulatory program");
  await expect(page.locator("#slide-59 aside.notes")).toHaveText("");
});
