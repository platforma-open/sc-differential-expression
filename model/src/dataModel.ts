import { kind } from "@platforma-open/milaboratories.sc-differential-expression.kind";
import { createPlDataTableStateV2, DataModelBuilder } from "@platforma-sdk/model";
import type { BlockArgs, BlockData, LegacyUiState } from "./types";

export const blockDataModel = new DataModelBuilder({ kind })
  .from<BlockData>("v1")
  // V1 kept these in two buckets and V3 keeps one. `anchorColumn` is dropped:
  // nothing in the UI ever wrote it.
  .upgradeLegacy<BlockArgs, LegacyUiState>(({ args, uiState }) => ({
    ...args,
    graphState: uiState.graphState,
    tableState: uiState.tableState,
  }))
  .init(({ params }) => ({
    countsRef: params?.countsRef,
    covariateRefs: params?.covariateRefs ?? [],
    contrastFactor: params?.contrastFactor,
    numerators: params?.numerators ?? [],
    denominator: params?.denominator,
    log2FCThreshold: params?.log2FCThreshold ?? 1,
    pAdjThreshold: params?.pAdjThreshold ?? 0.05,

    graphState: {
      title: "Differential gene expression",
      template: "dots",
      currentTab: null,
    },
    tableState: createPlDataTableStateV2(),
  }));
