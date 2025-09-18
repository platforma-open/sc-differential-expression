<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
import type { PlRef } from '@platforma-sdk/model';
import { plRefsEqual } from '@platforma-sdk/model';
import { PlAccordionSection, PlAgDataTableV2, PlAlert, PlBlockPage, PlBtnGhost, PlDropdown, PlDropdownMulti, PlDropdownRef, PlMaskIcon24, PlNumberField, PlRow, PlSlideModal, usePlDataTableSettingsV2 } from '@platforma-sdk/ui-vue';
import { computed, reactive } from 'vue';
import { useApp } from '../app';

const app = useApp();

const tableSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.topTableFilteredPt,
});

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
    !app.model.args.numerators.includes(op.value));
});

</script>

<template>
  <PlBlockPage>
    <template #title>Single-Cell Differential Expression Analysis</template>
    <template #append>
      <PlBtnGhost @click.stop="() => data.settingsOpen = true">
        Settings
        <template #append>
          <PlMaskIcon24 name="settings" />
        </template>
      </PlBtnGhost>
    </template>
    <PlAgDataTableV2
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
        required
        clearable @update:model-value="setInput"
      />
      <PlDropdownMulti v-model="app.model.args.covariateRefs" :options="covariateOptions" label="Covariates" required>
        <template #tooltip>
          <div>
            <strong>Covariates for differential expression</strong><br/>
            Select metadata variables to include in the differential expression model. Covariates help control for confounding factors (e.g., batch effects, cell cycle, technical variables) and improve the accuracy of differential expression results.
          </div>
        </template>
      </PlDropdownMulti>
      <PlDropdown v-model="app.model.args.contrastFactor" :options="contrastFactorOptions" label="Contrast factor" required>
        <template #tooltip>
          <div>
            <strong>Primary factor for comparison</strong><br/>
            Choose the main categorical variable for differential expression analysis. This defines the groups you want to compare (e.g., cell types, treatment conditions, disease states).
          </div>
        </template>
      </PlDropdown>
      <PlDropdownMulti v-model="app.model.args.numerators" :options="numeratorOptions" label="Numerator (test group)" required>
        <template #tooltip>
          <div>
            <strong>Test condition or group</strong><br/>
            Select the experimental condition or group of interest. Genes upregulated in this group will have positive log2 fold-changes in the results.
          </div>
        </template>
      </PlDropdownMulti>
      <PlDropdown v-model="app.model.args.denominator" :options="denominatorOptions" label="Denominator (reference group)" required>
        <template #tooltip>
          <div>
            <strong>Reference or control group</strong><br/>
            Select the baseline condition or control group for comparison. Expression levels in the numerator group will be compared relative to this reference group.
          </div>
        </template>
      </PlDropdown>

      <PlAccordionSection label="SIGNIFICANCE THRESHOLDS">
        <PlRow>
          <PlNumberField
            v-model="app.model.args.log2FCThreshold"
            label="Log2(FC)" :minValue="0" :step="0.1"
          >
            <template #tooltip>
              <div>
                <strong>Log2 fold-change threshold</strong><br/>
                Minimum absolute log2 fold-change required for a gene to be considered differentially expressed. Higher values (e.g., 1.0-2.0) identify more dramatically changed genes, while lower values (e.g., 0.5-0.8) capture subtler expression differences.
              </div>
            </template>
          </PlNumberField>
          <PlNumberField
            v-model="app.model.args.pAdjThreshold"
            label="Adjusted p-value" :minValue="0" :maxValue="1" :step="0.01"
          >
            <template #tooltip>
              <div>
                <strong>Statistical significance threshold</strong><br/>
                Maximum adjusted p-value for differential expression significance. Standard values are 0.05 (5% false discovery rate) or 0.01 (1% false discovery rate) for more stringent filtering.
              </div>
            </template>
          </PlNumberField>
        </PlRow>
        <!-- Add warnings if selected threshold are out of most commonly used bounds -->
        <PlAlert v-if="app.model.args.pAdjThreshold > 0.05" type="warn">
          {{ "Warning: The selected adjusted p-value threshold is higher than the commonly recommended 0.05" }}
        </PlAlert>
        <PlAlert v-if="app.model.args.log2FCThreshold < 0.6" type="warn">
          {{ "Warning: The selected Log2(FC) threshold may be too low for identifying robust differentially expressed genes" }}
        </PlAlert>
      </PlAccordionSection>
    </PlSlideModal>
  </PlBlockPage>
</template>
