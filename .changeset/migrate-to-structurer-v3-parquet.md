---
'@platforma-open/milaboratories.sc-differential-expression.model': minor
'@platforma-open/milaboratories.sc-differential-expression.ui': minor
'@platforma-open/milaboratories.sc-differential-expression': minor
'@platforma-open/milaboratories.sc-differential-expression.workflow': patch
'@platforma-open/milaboratories.sc-differential-expression.software': patch
---

Migrate the block onto the structurer (block-tools 2.16.12) with a full SDK upgrade: model and ui-vue 1.84.17, workflow-tengo 6.15.0, tengo-builder 4.1.11, package-builder 3.16.0. Adopts the canonical tool-managed layout, the slim facade, and the mandatory `kind/` package.

Move the model to BlockModelV3. SDK 1.83 removed the V1 `BlockModel`. Existing projects are upgraded once by `.upgradeLegacy`, and the args keep the V1 keys, so opening a project does not make it stale. The kind's init-params contract is every user setting (dataset, covariates, contrast factor, numerators, denominator, both thresholds), all optional. The UI no longer clears the chosen groups from a watcher on fetched values: changing the contrast factor clears them, and a group that disappears upstream shows a warning.

Exchange the counts matrix and the DE result tables as Parquet instead of CSV, and give the counts export and the Wilcoxon step 32 GiB / 4 CPU. The software builds the count matrix from categorical codes instead of per-row strings, the same memory fix as in cell-type-annotation. Results are unchanged.
