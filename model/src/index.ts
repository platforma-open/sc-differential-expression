import { kind } from "@platforma-open/milaboratories.sc-differential-expression.kind";
import type { InferOutputsType, PColumnIdAndSpec, PFrameHandle } from "@platforma-sdk/model";
import {
  BlockModelV3,
  createPFrameForGraphs,
  createPlDataTableV2,
  isPColumnSpec,
} from "@platforma-sdk/model";
import { blockDataModel } from "./dataModel";
import type { BlockArgs } from "./types";

export { blockDataModel } from "./dataModel";
export type { BlockArgs, BlockData } from "./types";

export const platforma = BlockModelV3.create({ dataModel: blockDataModel, kind })

  // The run gate. Throwing shows the reason in the UI. Replaces V1's
  // `.argsValid`. The projection strips the view state and returns the same
  // keys V1 stored as args.
  .args<BlockArgs>((data): BlockArgs => {
    const { graphState: _graphState, tableState: _tableState, ...args } = data;

    if (args.countsRef === undefined) {
      throw new Error("Dataset is required");
    }
    if (args.contrastFactor === undefined) {
      throw new Error("Contrast factor is required");
    }
    if (args.numerators.length === 0) {
      throw new Error("At least one numerator group is required");
    }
    if (args.denominator === undefined) {
      throw new Error("Denominator group is required");
    }

    return args;
  })

  // Inverse of the kind's init-params contract. `title` is left out: the UI
  // derives it from the dataset label.
  .templateParams((data) => ({
    countsRef: data.countsRef,
    covariateRefs: data.covariateRefs,
    contrastFactor: data.contrastFactor,
    numerators: data.numerators,
    denominator: data.denominator,
    log2FCThreshold: data.log2FCThreshold,
    pAdjThreshold: data.pAdjThreshold,
  }))

  .output("countsOptions", (ctx) =>
    ctx.resultPool.getOptions(
      (spec) =>
        isPColumnSpec(spec) &&
        spec.name === "pl7.app/rna-seq/countMatrix" &&
        spec.domain?.["pl7.app/rna-seq/normalized"] === "false",
      { includeNativeLabel: false, addLabelAsSuffix: true },
    ),
  )

  .output("metadataOptions", (ctx) =>
    ctx.resultPool.getOptions(
      (spec) =>
        isPColumnSpec(spec) &&
        (spec.name === "pl7.app/metadata" ||
          // @TODO: We will have a specific annotation in the future
          spec.name === "pl7.app/rna-seq/leidencluster" ||
          spec.name === "pl7.app/rna-seq/cellType"),
    ),
  )

  .output("numeratorOptions", (ctx) => {
    if (!ctx.data.contrastFactor) return undefined;

    const pColumn = ctx.resultPool.getPColumnByRef(ctx.data.contrastFactor);
    if (!pColumn) return undefined;

    return ctx.createPFrame([pColumn]);
  })

  // `withStatus` on the table and the volcano frame: PlAgDataTableV2 (via
  // usePlDataTableSettingsV2) and GraphMaker both take the status-wrapped
  // output and render the pending and error states from it.
  .outputWithStatus("topTableFilteredPt", (ctx) => {
    const pCols = ctx.outputs?.resolve("topTableFilteredPf")?.getPColumns();
    if (pCols === undefined) {
      return undefined;
    }

    return createPlDataTableV2(ctx, pCols, ctx.data.tableState);
  })

  .output("deErrors", (ctx) => {
    const pCols = ctx.outputs?.resolve("deErrors")?.getPColumns();
    if (pCols === undefined) {
      return undefined;
    }

    return ctx.createPFrame(pCols);
  })

  .outputWithStatus("topTablePf", (ctx): PFrameHandle | undefined => {
    const pCols = ctx.outputs?.resolve("topTablePf")?.getPColumns();
    if (pCols === undefined) {
      return undefined;
    }
    return createPFrameForGraphs(ctx, pCols);
  })

  .output("topTablePcols", (ctx) => {
    const pCols = ctx.outputs?.resolve("topTablePf")?.getPColumns();

    if (pCols === undefined || pCols.length === 0) {
      return undefined;
    }

    return pCols.map(
      (c) =>
        ({
          columnId: c.id,
          spec: c.spec,
        }) satisfies PColumnIdAndSpec,
    );
  })

  .output("isRunning", (ctx) => ctx.outputs?.getIsReadyOrError() === false)

  .sections((_ctx) => [
    { type: "link", href: "/", label: "Main" },
    { type: "link", href: "/volcano", label: "Volcano plot" },
  ])

  .title((ctx) =>
    ctx.data.title ? `Differential Expression - ${ctx.data.title}` : "Differential Expression",
  )

  .done();

export type BlockOutputs = InferOutputsType<typeof platforma>;
