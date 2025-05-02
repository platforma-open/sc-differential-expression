<script setup lang="ts">
import type { GraphMakerProps } from '@milaboratories/graph-maker';
import { GraphMaker } from '@milaboratories/graph-maker';
import '@milaboratories/graph-maker/styles';
import { useApp } from '../app';

const app = useApp();

const defaultOptions: GraphMakerProps['defaultOptions'] = [
  {
    inputName: 'x',
    selectedSource: {
      kind: 'PColumn',
      name: 'pl7.app/rna-seq/log2foldchange',
      valueType: 'Double',
      axesSpec: [
        {
          name: 'pl7.app/rna-seq/contrastGroup',
          type: 'String',
        },
        {
          name: 'pl7.app/rna-seq/geneId',
          type: 'String',
        },
      ],
    },
  },
  {
    inputName: 'y',
    selectedSource: {
      kind: 'PColumn',
      name: 'pl7.app/rna-seq/minlog10padj',
      valueType: 'Double',
      axesSpec: [
        {
          name: 'pl7.app/rna-seq/contrastGroup',
          type: 'String',
        },
        {
          name: 'pl7.app/rna-seq/geneId',
          type: 'String',
        },
      ],
    },
  },
  {
    inputName: 'grouping',
    selectedSource: {
      kind: 'PColumn',
      name: 'pl7.app/rna-seq/regulationDirection',
      valueType: 'String',
      axesSpec: [
        {
          name: 'pl7.app/rna-seq/contrastGroup',
          type: 'String',
        },
        {
          name: 'pl7.app/rna-seq/geneId',
          type: 'String',
        },
      ],
    },
  },
  {
    inputName: 'label',
    selectedSource: {
      name: 'pl7.app/rna-seq/geneId',
      type: 'String',
    },
  },
  {
    inputName: 'tooltipContent',
    selectedSource: {
      name: 'pl7.app/rna-seq/geneId',
      type: 'String',
    },
  },
];

</script>

<template>
  <GraphMaker
    v-model="app.model.ui.graphState" chartType="scatterplot-umap"
    :p-frame="app.model.outputs.topTablePf" :defaultOptions="defaultOptions"
  />
</template>
