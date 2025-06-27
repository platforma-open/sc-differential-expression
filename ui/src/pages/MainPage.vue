<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
import { PlAgDataTable, PlAgDataTableToolsPanel, PlNumberField, PlBlockPage, PlBtnGhost, PlDropdownRef, PlMaskIcon24, PlSlideModal, PlRow, PlAlert, PlDropdownMulti, PlDropdown, PlAccordionSection } from '@platforma-sdk/ui-vue';
import type { PlDataTableSettings } from '@platforma-sdk/ui-vue';
import { useApp } from '../app';
import { computed, reactive } from 'vue';
import type { PlRef } from '@platforma-sdk/model';
import { plRefsEqual } from '@platforma-sdk/model';

const app = useApp();

const tableSettings = computed<PlDataTableSettings>(() => ({
  sourceType: 'ptable',
  pTable: app.model.outputs.topTableFilteredPt,
}));

const data = reactive<{
  settingsOpen: boolean;
}>({
  settingsOpen: app.model.args.countsRef === undefined,
});

function setInput(inputRef?: PlRef) {
  app.model.args.countsRef = inputRef;
  if (inputRef)
    app.model.args.title = app.model.outputs.countsOptions?.find((o) => plRefsEqual(o.ref, inputRef))?.label;
  else
    app.model.args.title = undefined;
}

const covariateOptions = computed(() => {
  return app.model.outputs.metadataOptions?.map((v) => ({
    value: v.ref,
    label: v.label,
  })) ?? [];
});

const contrastFactorOptions = computed(() => {
  return app.model.args.covariateRefs.map((ref) => ({
    value: ref,
    label: covariateOptions.value.find((m) => m.value.name === ref.name)?.label ?? '',
  }));
});

const numeratorOptions = computed(() => {
  return app.model.outputs.denominatorOptions?.map((v) => ({
    value: v,
    label: v,
  }));
});

const denominatorOptions = computed(() => {
  return numeratorOptions.value?.filter((op) =>
    !app.model.args.numerator?.includes(op.value));
});

</script>

<template>
  <PlBlockPage>
    <template #title>sc Differential Expression</template>
    <template #append>
      <!-- PlAgDataTableToolsPanel controls showing  Export column and filter-->
      <PlAgDataTableToolsPanel/>
      <PlBtnGhost @click.stop="() => data.settingsOpen = true">
        Settings
        <template #append>
          <PlMaskIcon24 name="settings" />
        </template>
      </PlBtnGhost>
    </template>
    <PlAgDataTable
      v-model="app.model.ui.tableState"
      :settings="tableSettings"
      show-columns-panel
      show-export-button
    />
    <PlSlideModal v-model="data.settingsOpen">
      <template #title>Settings</template>
      <PlDropdownRef
        v-model="app.model.args.countsRef" :options="app.model.outputs.countsOptions"
        label="Select dataset"
        clearable @update:model-value="setInput"
      />
      <PlDropdownMulti v-model="app.model.args.covariateRefs" :options="covariateOptions" label="Covariates" />
      <PlDropdown v-model="app.model.args.contrastFactor" :options="contrastFactorOptions" label="Contrast factor" />
      <PlDropdown v-model="app.model.args.numerator" :options="numeratorOptions" label="Numerator" />
      <PlDropdown v-model="app.model.args.denominator" :options="denominatorOptions" label="Denominator" />

      <PlAccordionSection label="THRESHOLD PARAMETERS">
        <PlRow>
          <PlNumberField
            v-model="app.model.args.log2FCThreshold"
            label="Log2(FC)" :minValue="0" :step="0.1"
          >
            <template #tooltip>
              Select a valid absolute log2(FC) threshold for identifying
              significant genes.
            </template>
          </PlNumberField>
          <PlNumberField
            v-model="app.model.args.pAdjThreshold"
            label="Adjusted p-value" :minValue="0" :maxValue="1" :step="0.01"
          />
        </PlRow>
        <!-- Add warnings if selected threshold are out of most commonly used bounds -->
        <PlAlert v-if="app.model.args.pAdjThreshold > 0.05" type="warn">
          {{ "Warning: The selected adjusted p-value threshold is higher than the most commonly used 0.05" }}
        </PlAlert>
        <PlAlert v-if="app.model.args.log2FCThreshold < 0.6" type="warn">
          {{ "Warning: The selected Log2(FC) threshold may be too low for most use cases" }}
        </PlAlert>
      </PlAccordionSection>
    </PlSlideModal>
  </PlBlockPage>
</template>
