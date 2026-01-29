# R Package Logging And Progress (Best Practice)

Use this guide when adding or modifying user-facing output, logging, or progress
reporting in this repository.

## Principles
- Keep console output concise and user-controlled.
- Never print from worker processes.
- Prefer structured, consistent messaging utilities.
- Make progress reporting robust in both single-core and multi-core runs.

## Messaging
- Prefer the shared timestamped logger in `R/utils_logging.R`.
- Use `.log_inform()` for informational messages.
- Use `.log_warn()` for warnings and `.log_abort()` for errors.
- Gate all informational output behind a `verbose` (or `quiet`) argument.
- Avoid `cat()` or `print()` for user-facing output.
- Use `pkg::fun` namespacing; avoid `library()`/`require()` in package code.

## Progress
- Single-core: use `cli` progress (e.g., `cli::cli_progress_*`) or simple
  `.log_inform()` checkpoints for coarse stages.
- Multi-core:
  - Do not emit progress/messages inside workers.
  - Prefer `progressr` with a parallel backend that supports it (e.g., `future`)
    so progress is aggregated in the main process.
  - If using `mclapply`/`parallel`, report progress per chunk only in the main
    process or write logs to files.

## Logging To Files
- When detailed progress is needed for multi-core tasks, add a `log_file`
  argument and write structured lines (e.g., TSV: timestamp, stage, item, status).
- Logging must not change outputs or performance-critical paths by default.

## Defaults
- Default: `verbose = TRUE` with concise `cli::cli_inform()` messages.
- Provide `verbose = FALSE` or `quiet = TRUE` to silence output.
- For long-running steps, emit a start and finish message; avoid per-item output
  in the console.

## Package Quality
- Add roxygen2 docs for exported functions and keep examples fast.
- Validate inputs with clear, actionable error messages.
- Keep R CMD check-friendly behavior: no hidden globals, no side effects on load,
  deterministic where needed (e.g., set a seed in stochastic examples).
- Avoid unnecessary complexity; document if a change increases complexity.
- Prefer clarity over cleverness; follow repository style.

## Plotting
- Prefer ggplot2.
- Use bold x/y axis titles.
- Include an informative title (what + condition + units).
- Add a caption when assumptions or filters are applied.
- If the package defines a theme helper (e.g., `theme_pkg()`), use it.

## Outputs
- When asked to return file contents, output complete files in fenced code blocks,
  and provide tests in a separate block if requested.

## Edit Boundaries
- If an edit-only mode is specified, modify code only between explicit markers
  (e.g., `# BEGIN EDIT` / `# END EDIT`) and keep other sections byte-identical.
