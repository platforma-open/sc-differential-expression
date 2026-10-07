import type { GraphMakerState } from "@milaboratories/graph-maker";
import type { PlDataTableStateV2, PlRef } from "@platforma-sdk/model";

/**
 * Workflow-facing args — every user-committed analysis decision. Unchanged
 * from the V1 shape on purpose: keeping these bytes identical is what stops
 * every existing project going stale the moment it is opened on this version.
 *
 * `title` is not read by the workflow. The UI sets it from the dataset label,
 * and `.title()` reads it.
 */
export type BlockArgs = {
  countsRef?: PlRef;
  covariateRefs: PlRef[];
  contrastFactor?: PlRef;
  denominator?: string;
  numerators: string[];
  log2FCThreshold: number;
  pAdjThreshold: number;
  title?: string;
};

/** The V1 view state, as it is on disk. Read only by `.upgradeLegacy`. */
export type LegacyUiState = {
  anchorColumn?: PlRef;
  graphState: GraphMakerState;
  tableState: PlDataTableStateV2;
};

/**
 * Unified V3 data — the args above plus the table and plot view state.
 *
 * Every settings-panel field is an analysis decision the user commits and then
 * presses Run for, so none of them belongs in `prerunArgs`. No V1 field was
 * bent to dodge the staleness gate, so the V1 split is already the split V3
 * wants.
 */
export type BlockData = BlockArgs & {
  graphState: GraphMakerState;
  tableState: PlDataTableStateV2;
};
