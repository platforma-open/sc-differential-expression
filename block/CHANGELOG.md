# @platforma-open/milaboratories.sc-differential-expression

## 1.8.1

### Patch Changes

- 6509608: Release every package through npm trusted publishing (OIDC) so each package's trusted publisher is used before it expires.

## 1.8.0

### Minor Changes

- e473d86: Migrate the block onto the structurer (block-tools 2.16.12) with a full SDK upgrade: model and ui-vue 1.84.17, workflow-tengo 6.15.0, tengo-builder 4.1.11, package-builder 3.16.0. Adopts the canonical tool-managed layout, the slim facade, and the mandatory `kind/` package.

  Move the model to BlockModelV3. SDK 1.83 removed the V1 `BlockModel`. Existing projects are upgraded once by `.upgradeLegacy`, and the args keep the V1 keys, so opening a project does not make it stale. The kind's init-params contract is every user setting (dataset, covariates, contrast factor, numerators, denominator, both thresholds), all optional. The UI no longer clears the chosen groups from a watcher on fetched values: changing the contrast factor clears them, and a group that disappears upstream shows a warning.

  Exchange the counts matrix and the DE result tables as Parquet instead of CSV, and give the counts export and the Wilcoxon step 32 GiB / 4 CPU. The software builds the count matrix from categorical codes instead of per-row strings, the same memory fix as in cell-type-annotation. Results are unchanged.

## 1.3.1

### Patch Changes

- @platforma-open/milaboratories.sc-differential-expression.workflow@1.7.1

## 1.3.0

### Minor Changes

- 9183a8c: Enable deduplication and block metadata update.

### Patch Changes

- Updated dependencies [9183a8c]
  - @platforma-open/milaboratories.sc-differential-expression.workflow@1.7.0

## 1.2.15

### Patch Changes

- Updated dependencies [dfcff01]
  - @platforma-open/milaboratories.sc-differential-expression.workflow@1.6.3

## 1.2.14

### Patch Changes

- @platforma-open/milaboratories.sc-differential-expression.workflow@1.6.2

## 1.2.13

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
  - @platforma-open/milaboratories.sc-differential-expression.ui@1.4.1
  - @platforma-open/milaboratories.sc-differential-expression.workflow@1.6.1

## 1.2.12

### Patch Changes

- Updated dependencies [498627e]
  - @platforma-open/milaboratories.sc-differential-expression.workflow@1.6.0
  - @platforma-open/milaboratories.sc-differential-expression.model@1.4.0
  - @platforma-open/milaboratories.sc-differential-expression.ui@1.4.0

## 1.2.11

### Patch Changes

- Updated dependencies [a2be6d6]
  - @platforma-open/milaboratories.sc-differential-expression.model@1.3.4
  - @platforma-open/milaboratories.sc-differential-expression.ui@1.3.4

## 1.2.10

### Patch Changes

- Updated dependencies [d560773]
  - @platforma-open/milaboratories.sc-differential-expression.workflow@1.5.4
  - @platforma-open/milaboratories.sc-differential-expression.model@1.3.3
  - @platforma-open/milaboratories.sc-differential-expression.ui@1.3.3

## 1.2.9

### Patch Changes

- Updated dependencies [ab9a5ce]
  - @platforma-open/milaboratories.sc-differential-expression.model@1.3.2
  - @platforma-open/milaboratories.sc-differential-expression.workflow@1.5.3
  - @platforma-open/milaboratories.sc-differential-expression.ui@1.3.2

## 1.2.8

### Patch Changes

- Updated dependencies [bf2d04f]
  - @platforma-open/milaboratories.sc-differential-expression.workflow@1.5.2

## 1.2.7

### Patch Changes

- d16afba: Fix merge publication

## 1.2.6

### Patch Changes

- Updated dependencies [18c0d47]
  - @platforma-open/milaboratories.sc-differential-expression.model@1.3.1
  - @platforma-open/milaboratories.sc-differential-expression.workflow@1.5.1
  - @platforma-open/milaboratories.sc-differential-expression.ui@1.3.1

## 1.2.5

### Patch Changes

- Updated dependencies [177cc06]
  - @platforma-open/milaboratories.sc-differential-expression.workflow@1.5.0

## 1.2.4

### Patch Changes

- @platforma-open/milaboratories.sc-differential-expression.workflow@1.4.1

## 1.2.3

### Patch Changes

- Updated dependencies [cfd865c]
  - @platforma-open/milaboratories.sc-differential-expression.workflow@1.4.0

## 1.2.2

### Patch Changes

- Updated dependencies [0f60e23]
  - @platforma-open/milaboratories.sc-differential-expression.workflow@1.3.0
  - @platforma-open/milaboratories.sc-differential-expression.model@1.3.0
  - @platforma-open/milaboratories.sc-differential-expression.ui@1.3.0

## 1.2.1

### Patch Changes

- Updated dependencies [f2c5fd6]
  - @platforma-open/milaboratories.sc-differential-expression.workflow@1.2.0
  - @platforma-open/milaboratories.sc-differential-expression.model@1.2.0
  - @platforma-open/milaboratories.sc-differential-expression.ui@1.2.0

## 1.2.0

### Minor Changes

- 5e8aecf: allow prepare venv on Windows

## 1.1.1

### Patch Changes

- Updated dependencies [4ffab78]
  - @platforma-open/milaboratories.sc-differential-expression.workflow@1.1.1

## 1.1.0

### Minor Changes

- eaf2947: Block first version

### Patch Changes

- Updated dependencies [eaf2947]
  - @platforma-open/milaboratories.sc-differential-expression.workflow@1.1.0
  - @platforma-open/milaboratories.sc-differential-expression.model@1.1.0
  - @platforma-open/milaboratories.sc-differential-expression.ui@1.1.0
