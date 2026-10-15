import { describe, expect, it } from "vitest";
import { kind } from "./index";

const parse = (params: unknown) => kind.parseInitializationParams(params);

/** Refs exactly as `resultPool.getOptions` mints them for the dropdowns. */
const COUNTS_REF = {
  __isRef: true,
  blockId: "cell-ranger-1",
  name: "pl7.app/rna-seq/countMatrix",
} as const;
const CONDITION_REF = {
  __isRef: true,
  blockId: "samples-and-data-1",
  name: "pl7.app/metadata",
} as const;

describe("the params envelope", () => {
  it.each([
    ["null", null],
    ["a number", 5],
    ["an array", ["numerators"]],
    ["a string", "denominator=control"],
  ])("rejects %s", (_label, value) => {
    expect(() => parse(value)).toThrow();
  });

  it("accepts an empty object — a block created with nothing seeded", () => {
    expect(parse({})).toEqual({});
  });

  it("drops keys the contract does not name", () => {
    expect(parse({ denominator: "control", title: "Differential Expression - X" })).toEqual({
      denominator: "control",
    });
  });

  it("returns a full set of params unchanged", () => {
    const params = {
      countsRef: COUNTS_REF,
      covariateRefs: [CONDITION_REF],
      contrastFactor: CONDITION_REF,
      numerators: ["treated"],
      denominator: "control",
      log2FCThreshold: 1,
      pAdjThreshold: 0.05,
    };
    expect(parse(params)).toEqual(params);
  });
});

describe("refs", () => {
  it.each([
    ["an object missing the ref marker", { blockId: "b1", name: "pl7.app/metadata" }],
    ["a bare column name", "pl7.app/metadata"],
    ["null", null],
  ])("rejects %s", (_label, ref) => {
    expect(() => parse({ countsRef: ref })).toThrow("'countsRef' must be");
    expect(() => parse({ contrastFactor: ref })).toThrow("'contrastFactor' must be");
    expect(() => parse({ covariateRefs: [ref] })).toThrow("'covariateRefs' must be");
  });
});

describe("groups", () => {
  it("rejects a numerator that is not a string", () => {
    expect(() => parse({ numerators: ["treated", 2] })).toThrow("'numerators' must be");
  });

  it("rejects a denominator that is not a string", () => {
    expect(() => parse({ denominator: ["control"] })).toThrow("'denominator' must be");
  });
});

describe("thresholds", () => {
  it.each([-0.1, Number.NaN, "1"])("rejects log2FCThreshold %s", (value) => {
    expect(() => parse({ log2FCThreshold: value })).toThrow("'log2FCThreshold' must be");
  });

  it.each([-0.01, 1.5, Number.POSITIVE_INFINITY])("rejects pAdjThreshold %s", (value) => {
    expect(() => parse({ pAdjThreshold: value })).toThrow("'pAdjThreshold' must be");
  });

  it("accepts the boundary values", () => {
    expect(parse({ log2FCThreshold: 0, pAdjThreshold: 1 })).toEqual({
      log2FCThreshold: 0,
      pAdjThreshold: 1,
    });
  });
});
