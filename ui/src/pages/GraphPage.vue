<script setup lang="ts">
import type { PredefinedGraphOption } from '@milaboratories/graph-maker';
import { GraphMaker } from '@milaboratories/graph-maker';
import '@milaboratories/graph-maker/styles';
import type { PColumnIdAndSpec } from '@platforma-sdk/model';
import { computed } from 'vue';
import { useApp } from '../app';

const app = useApp();

// const defaultOptions = computed((): GraphMakerProps['defaultOptions'] => {
const defaultOptions = computed((): PredefinedGraphOption<'scatterplot-umap'>[] | undefined => {
  if (!app.model.outputs.topTablePcols)
    return undefined;

  const topTablePcols = app.model.outputs.topTablePcols;
  function getIndex(name: string, pcols: PColumnIdAndSpec[]): number {
    return pcols.findIndex((p) => (p.spec.name === name));
  }

  // Check if all required columns exist
  const requiredColumns = [
    'pl7.app/rna-seq/log2foldchange',
    'pl7.app/rna-seq/minlog10padj',
    'pl7.app/rna-seq/regulationDirection',
  ];

  for (const columnName of requiredColumns) {
    if (getIndex(columnName, topTablePcols) === -1) {
      return undefined;
    }
  }

  const defaults: PredefinedGraphOption<'scatterplot-umap'>[] = [
    {
      inputName: 'x',
      selectedSource: topTablePcols[getIndex('pl7.app/rna-seq/log2foldchange', topTablePcols)].spec,
    },
    {
      inputName: 'y',
      selectedSource: topTablePcols[getIndex('pl7.app/rna-seq/minlog10padj', topTablePcols)].spec,
    },
    {
      inputName: 'grouping',
      selectedSource: topTablePcols[getIndex('pl7.app/rna-seq/regulationDirection', topTablePcols)].spec,
    },
    {
      inputName: 'label',
      selectedSource: topTablePcols[getIndex('pl7.app/rna-seq/regulationDirection', topTablePcols)].spec.axesSpec[1],
    },
    {
      inputName: 'tooltipContent',
      selectedSource: topTablePcols[getIndex('pl7.app/rna-seq/regulationDirection', topTablePcols)].spec.axesSpec[1],
    },
  ];
  return defaults;
});

</script>

<template>
  <GraphMaker
    v-model="app.model.ui.graphState" chartType="scatterplot-umap"
    :p-frame="app.model.outputs.topTablePf" :defaultOptions="defaultOptions"
  />
</template>
