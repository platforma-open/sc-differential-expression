# @platforma-open/milaboratories.sc-differential-expression.workflow

## 1.8.0

### Patch Changes

- e473d86: Migrate the block onto the structurer (block-tools 2.16.12) with a full SDK upgrade: model and ui-vue 1.84.17, workflow-tengo 6.15.0, tengo-builder 4.1.11, package-builder 3.16.0. Adopts the canonical tool-managed layout, the slim facade, and the mandatory `kind/` package.

  Move the model to BlockModelV3. SDK 1.83 removed the V1 `BlockModel`. Existing projects are upgraded once by `.upgradeLegacy`, and the args keep the V1 keys, so opening a project does not make it stale. The kind's init-params contract is every user setting (dataset, covariates, contrast factor, numerators, denominator, both thresholds), all optional. The UI no longer clears the chosen groups from a watcher on fetched values: changing the contrast factor clears them, and a group that disappears upstream shows a warning.

  Exchange the counts matrix and the DE result tables as Parquet instead of CSV, and give the counts export and the Wilcoxon step 32 GiB / 4 CPU. The software builds the count matrix from categorical codes instead of per-row strings, the same memory fix as in cell-type-annotation. Results are unchanged.

- Updated dependencies [e473d86]
  - @platforma-open/milaboratories.sc-differential-expression.software@1.3.4

## 1.7.1

### Patch Changes

- Updated dependencies [e020722]
  - @platforma-open/milaboratories.sc-differential-expression.software@1.3.3

## 1.7.0

### Minor Changes

- 9183a8c: Enable deduplication and block metadata update.

## 1.6.3

### Patch Changes

- dfcff01: Update to Parquet

## 1.6.2

### Patch Changes

- Updated dependencies [41bdf7b]
  - @platforma-open/milaboratories.sc-differential-expression.software@1.3.2

## 1.6.1

### Patch Changes

- 9a84bf7: technical release
- 89faa47: technical release
- 03d124f: technical release
- e8e3c63: technical release
- Updated dependencies [9a84bf7]
- Updated dependencies [89faa47]
- Updated dependencies [03d124f]
- Updated dependencies [e8e3c63]
  - @platforma-open/milaboratories.sc-differential-expression.software@1.3.1

## 1.6.0

### Minor Changes

- 498627e: Allow cell annotations in contrast and show warnings if there are issues with the analysis

### Patch Changes

- Updated dependencies [498627e]
  - @platforma-open/milaboratories.sc-differential-expression.software@1.3.0

## 1.5.4

### Patch Changes

- d560773: Refactoring and allow multiple numerators

## 1.5.3

### Patch Changes

- Updated dependencies [ab9a5ce]
  - @platforma-open/milaboratories.sc-differential-expression.software@1.2.2

## 1.5.2

### Patch Changes

- bf2d04f: Updated trace

## 1.5.1

### Patch Changes

- Updated dependencies [18c0d47]
  - @platforma-open/milaboratories.sc-differential-expression.software@1.2.1

## 1.5.0

### Minor Changes

- 177cc06: Add batch support

## 1.4.1

### Patch Changes

- Updated dependencies [5889cc3]
  - @platforma-open/milaboratories.sc-differential-expression.software@1.2.0

## 1.4.0

### Minor Changes

- cfd865c: Add trace to exports

## 1.3.0

### Minor Changes

- 0f60e23: Update SDK and migrate to PlAgDataTableV2

## 1.2.0

### Minor Changes

- f2c5fd6: Added gene symbols for better interpretation of the results.

## 1.1.1

### Patch Changes

- 4ffab78: chore: update deps

## 1.1.0

### Minor Changes

- eaf2947: Block first version

### Patch Changes

- Updated dependencies [eaf2947]
  - @platforma-open/milaboratories.sc-differential-expression.software@1.1.0
