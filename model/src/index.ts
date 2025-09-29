import type { GraphMakerState } from '@milaboratories/graph-maker';
import type {
  InferOutputsType,
  PColumnIdAndSpec,
  PFrameHandle,
  PlDataTableStateV2,
  PlRef,
} from '@platforma-sdk/model';
import {
  BlockModel,
  createPFrameForGraphs,
  createPlDataTableV2,
  isPColumnSpec,
} from '@platforma-sdk/model';

export type UiState = {
  anchorColumn?: PlRef;
  graphState: GraphMakerState;
  tableState: PlDataTableStateV2;
};

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

export const model = BlockModel.create()

  .withArgs<BlockArgs>({
    covariateRefs: [],
    numerators: [],
    log2FCThreshold: 1,
    pAdjThreshold: 0.05,
  })

  .withUiState<UiState>({
    tableState: {
      gridState: {},
      pTableParams: {
        sorting: [],
        filters: [],
      },
    },
    graphState: {
      title: 'Differential gene expression',
      template: 'dots',
      currentTab: null,
    },
  })

  .argsValid((ctx) => (
    (ctx.args.numerators !== undefined)) && (ctx.args.numerators.length !== 0)
  && (ctx.args.denominator !== undefined)
  && (ctx.args.log2FCThreshold !== undefined)
  && (ctx.args.pAdjThreshold !== undefined),
  )

  .output('countsOptions', (ctx) =>
    ctx.resultPool.getOptions((spec) => isPColumnSpec(spec)
      && spec.name === 'pl7.app/rna-seq/countMatrix'
      && spec.domain?.['pl7.app/rna-seq/normalized'] === 'false'
      // && spec.annotations?.['pl7.app/hideDataFromGraphs'] === 'true'
    , { includeNativeLabel: false, addLabelAsSuffix: true }),
  )

  .output('metadataOptions', (ctx) =>
    ctx.resultPool.getOptions((spec) => isPColumnSpec(spec)
      && ((spec.name === 'pl7.app/metadata')
        // @TODO: We will have a specific annotation in the future
        || (spec.name === 'pl7.app/rna-seq/leidencluster')
        || (spec.name === 'pl7.app/rna-seq/cellType'))),
  )

  .output('numeratorOptions', (ctx) => {
    if (!ctx.args.contrastFactor) return undefined;

    const pColumn = ctx.resultPool.getPColumnByRef(ctx.args.contrastFactor);
    if (!pColumn) return undefined;

    return ctx.createPFrame([pColumn]);
  })

  .output('topTableFilteredPt', (ctx) => {
    const pCols = ctx.outputs?.resolve('topTableFilteredPf')?.getPColumns();
    if (pCols === undefined) {
      return undefined;
    }

    return createPlDataTableV2(ctx, pCols, ctx.uiState?.tableState);
  })

  .output('topTablePf', (ctx): PFrameHandle | undefined => {
    const pCols = ctx.outputs?.resolve('topTablePf')?.getPColumns();
    if (pCols === undefined) {
      return undefined;
    }
    return createPFrameForGraphs(ctx, pCols);
  })

  .output('topTablePcols', (ctx) => {
    const pCols = ctx.outputs?.resolve('topTablePf')?.getPColumns();

    if (pCols === undefined || pCols.length === 0) {
      return undefined;
    }

    return pCols.map(
      (c) =>
        ({
          columnId: c.id,
          spec: c.spec,
        } satisfies PColumnIdAndSpec),
    );
  })

  .output('anchorSpec', (ctx) => {
    // return the Reference of the p-column selected as input dataset in Settings
    if (!ctx.uiState?.anchorColumn) return undefined;

    // Get the specs of that selected p-column
    const anchorColumn = ctx.resultPool.getPColumnByRef(ctx.uiState?.anchorColumn);
    const anchorSpec = anchorColumn?.spec;
    if (!anchorSpec) {
      console.error('Anchor spec is undefined or is not PColumnSpec', anchorSpec);
      return undefined;
    }

    return anchorSpec;
  })

  .output('isRunning', (ctx) => ctx.outputs?.getIsReadyOrError() === false)

  .sections((_ctx) => ([
    { type: 'link', href: '/', label: 'Main' },
    { type: 'link', href: '/volcano', label: 'Volcano plot' },
  ]))

  .title((ctx) =>
    ctx.args.title
      ? `Differential Expression - ${ctx.args.title}`
      : 'Differential Expression',
  )

  .done(2);

export type BlockOutputs = InferOutputsType<typeof model>;
