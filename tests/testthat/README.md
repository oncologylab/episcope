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
  - `topic` for topic benchmark/report tests.
  - `benchmark` for benchmark-script coverage.
  - `tests` for tests that validate the test suite itself.
  - `utils`, `golem`, `fct`, or `app` for shared package helpers.

## Examples

- `test-module1-predict-tfbs.R`
- `test-module2-link-tf-targets.R`
- `test-module3-topic-models.R`
- `test-topic-method-k-html-reports.R`
- `test-tests-naming.R`
- `helper-module1.R`
