import { assertParamsObject, defineBlockKind } from "@platforma-sdk/block-kind";
import type { PlRef } from "@platforma-sdk/model";
import { isPlRef } from "@platforma-sdk/model";
import { name, version } from "../package.json" with { type: "json" };

/**
 * This block's init-params contract — the settings a user picks by hand: the
 * raw counts matrix, the covariates, the contrast factor with its numerator and
 * denominator groups, and the two significance thresholds.
 *
 * `title` is absent. The UI sets it from the label of the chosen dataset, so
 * there is nothing for a template to carry.
 *
 * Every field is optional. A block with no dataset chosen is ordinary state the
 * UI reaches — it is what a freshly created block looks like — and the
 * projection hands that state back untouched, so a required field would make
 * this block export a file its own kind refuses.
 */
export type BlockParams = {
  countsRef?: PlRef;
  covariateRefs?: PlRef[];
  contrastFactor?: PlRef;
  numerators?: string[];
  denominator?: string;
  log2FCThreshold?: number;
  pAdjThreshold?: number;
};

// Identity (`name`/`version`) comes from this package's own `package.json`, so
// the on-wire `{name}@{version}` reference can never drift from what npm
// publishes; the bundler inlines the JSON import.
export const kind = defineBlockKind<BlockParams>({
  name,
  version,
  parseInitializationParams,
});

// Internals

/** The same contract at runtime, for params arriving from a template file rather than typed code. */
function parseInitializationParams(value: unknown): BlockParams {
  assertParamsObject(value);

  const {
    countsRef,
    covariateRefs,
    contrastFactor,
    numerators,
    denominator,
    log2FCThreshold,
    pAdjThreshold,
  } = value;

  if (countsRef !== undefined && !isPlRef(countsRef)) {
    throw new Error(refMessage("countsRef"));
  }
  if (
    covariateRefs !== undefined &&
    !(Array.isArray(covariateRefs) && covariateRefs.every(isPlRef))
  ) {
    throw new Error("'covariateRefs' must be a list of references to upstream columns.");
  }
  if (contrastFactor !== undefined && !isPlRef(contrastFactor)) {
    throw new Error(refMessage("contrastFactor"));
  }
  if (numerators !== undefined && !isStringArray(numerators)) {
    throw new Error("'numerators' must be a list of group names.");
  }
  if (denominator !== undefined && typeof denominator !== "string") {
    throw new Error("'denominator' must be a group name.");
  }
  // A floor on |log2 fold change|, so any non-negative number.
  if (log2FCThreshold !== undefined && !isNonNegativeNumber(log2FCThreshold)) {
    throw new Error("'log2FCThreshold' must be a number greater than or equal to 0.");
  }
  // A ceiling on the adjusted p-value, so a probability.
  if (pAdjThreshold !== undefined && !isProbability(pAdjThreshold)) {
    throw new Error("'pAdjThreshold' must be a number between 0 and 1.");
  }

  return {
    countsRef,
    covariateRefs,
    contrastFactor,
    numerators,
    denominator,
    log2FCThreshold,
    pAdjThreshold,
  };
}

function refMessage(field: string): string {
  return `'${field}' must be a reference to an upstream column, written as { __isRef: true, blockId, name }.`;
}

function isStringArray(value: unknown): value is string[] {
  return Array.isArray(value) && value.every((v) => typeof v === "string");
}

function isNonNegativeNumber(value: unknown): value is number {
  return typeof value === "number" && Number.isFinite(value) && value >= 0;
}

function isProbability(value: unknown): value is number {
  return typeof value === "number" && Number.isFinite(value) && value >= 0 && value <= 1;
}
