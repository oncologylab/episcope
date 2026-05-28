# Test File Naming

Use kebab-case for all `testthat` files.

## Rules

- Test files must use `test-<area>-<feature>.R`.
- Helper files must use `helper-<area>.R`.
- Use lowercase words separated by hyphens.
- Avoid underscores in test file names.
- Keep the area prefix aligned with the package structure:
  - `module1` for TFBS prediction tests.
  - `module2` for TF-target linking and report tests.
  - `module3` for differential GRN and topic-model tests.
  - `topic` for topic-model benchmark and report tests.
  - `benchmark` for benchmark-script coverage.
  - `tests` for tests that validate the test suite itself.
  - `app` for Shiny application contract and UI/server helpers.
  - `config` for project configuration behavior.

## Examples

- `test-module1-predict-tfbs.R`
- `test-module1-footprints.R`
- `test-module2-link-tf-targets.R`
- `test-app-ui-utils.R`
- `test-module3-topic-models.R`
- `test-topic-reports.R`
- `test-test-suite-naming.R`
- `helper-module1.R`
