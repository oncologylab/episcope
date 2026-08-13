(() => {
  function closeDemo(demo) {
    if (!demo) return;
    const frame = demo.querySelector(".demo-frame");
    const activate = demo.querySelector(".demo-activate");
    const fallback = demo.querySelector(".demo-fallback");
    frame.hidden = true;
    activate.hidden = false;
    fallback.hidden = false;
    activate.focus({ preventScroll: true });
    if (window.Reveal) Reveal.focus();
  }

  document.querySelectorAll(".demo").forEach((demo) => {
    const activate = demo.querySelector(".demo-activate");
    const fallback = demo.querySelector(".demo-fallback");
    const frame = demo.querySelector(".demo-frame");
    const iframe = demo.querySelector("iframe");
    const status = demo.querySelector(".demo-status");
    activate.addEventListener("click", () => {
      activate.hidden = true;
      fallback.hidden = true;
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
