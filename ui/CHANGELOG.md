# @platforma-open/milaboratories.sc-differential-expression.ui

## 1.8.1

### Patch Changes

- 6509608: Release every package through npm trusted publishing (OIDC) so each package's trusted publisher is used before it expires.
- Updated dependencies [6509608]
  - @platforma-open/milaboratories.sc-differential-expression.model@1.8.1

## 1.8.0

### Minor Changes

- e473d86: Migrate the block onto the structurer (block-tools 2.16.12) with a full SDK upgrade: model and ui-vue 1.84.17, workflow-tengo 6.15.0, tengo-builder 4.1.11, package-builder 3.16.0. Adopts the canonical tool-managed layout, the slim facade, and the mandatory `kind/` package.

  Move the model to BlockModelV3. SDK 1.83 removed the V1 `BlockModel`. Existing projects are upgraded once by `.upgradeLegacy`, and the args keep the V1 keys, so opening a project does not make it stale. The kind's init-params contract is every user setting (dataset, covariates, contrast factor, numerators, denominator, both thresholds), all optional. The UI no longer clears the chosen groups from a watcher on fetched values: changing the contrast factor clears them, and a group that disappears upstream shows a warning.

  Exchange the counts matrix and the DE result tables as Parquet instead of CSV, and give the counts export and the Wilcoxon step 32 GiB / 4 CPU. The software builds the count matrix from categorical codes instead of per-row strings, the same memory fix as in cell-type-annotation. Results are unchanged.

### Patch Changes

- Updated dependencies [e473d86]
  - @platforma-open/milaboratories.sc-differential-expression.model@1.8.0

## 1.4.1

### Patch Changes

- 9a84bf7: technical release
- 89faa47: technical release
- 03d124f: technical release
- e8e3c63: technical release
- Updated dependencies [9a84bf7]
- Updated dependencies [89faa47]
- Updated dependencies [03d124f]
- Updated dependencies [e8e3c63]
  - @platforma-open/milaboratories.sc-differential-expression.model@1.4.1

## 1.4.0

### Minor Changes

- 498627e: Allow cell annotations in contrast and show warnings if there are issues with the analysis

### Patch Changes

- Updated dependencies [498627e]
  - @platforma-open/milaboratories.sc-differential-expression.model@1.4.0

## 1.3.4

### Patch Changes

- Updated dependencies [a2be6d6]
  - @platforma-open/milaboratories.sc-differential-expression.model@1.3.4

## 1.3.3

### Patch Changes

- d560773: Refactoring and allow multiple numerators
- Updated dependencies [d560773]
  - @platforma-open/milaboratories.sc-differential-expression.model@1.3.3

## 1.3.2

### Patch Changes

- Updated dependencies [ab9a5ce]
  - @platforma-open/milaboratories.sc-differential-expression.model@1.3.2

## 1.3.1

### Patch Changes

- Updated dependencies [18c0d47]
  - @platforma-open/milaboratories.sc-differential-expression.model@1.3.1

## 1.3.0

### Minor Changes

- 0f60e23: Update SDK and migrate to PlAgDataTableV2

### Patch Changes

- Updated dependencies [0f60e23]
  - @platforma-open/milaboratories.sc-differential-expression.model@1.3.0

## 1.2.0

### Minor Changes

- f2c5fd6: Added gene symbols for better interpretation of the results.

### Patch Changes

- Updated dependencies [f2c5fd6]
  - @platforma-open/milaboratories.sc-differential-expression.model@1.2.0

## 1.1.0

### Minor Changes

- eaf2947: Block first version

### Patch Changes

- Updated dependencies [eaf2947]
  - @platforma-open/milaboratories.sc-differential-expression.model@1.1.0
