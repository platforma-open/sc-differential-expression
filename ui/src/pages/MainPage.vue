<script setup lang="ts">
import type { PlRef } from "@platforma-sdk/model";
import { PFrameImpl, plRefsEqual } from "@platforma-sdk/model";
import {
  PlAccordionSection,
  PlAgDataTableV2,
  PlAlert,
  PlBlockPage,
  PlBtnGhost,
  PlDropdown,
  PlDropdownMulti,
  PlDropdownRef,
  PlMaskIcon24,
  PlNumberField,
  PlRow,
  PlSlideModal,
  usePlDataTableSettingsV2,
  useWatchFetch,
} from "@platforma-sdk/ui-vue";
import { computed, reactive } from "vue";
import { useApp } from "../app";

const app = useApp();

const tableSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.topTableFilteredPt,
});

const data = reactive<{
  settingsOpen: boolean;
}>({
  settingsOpen: app.model.data.countsRef === undefined,
});

function setInput(inputRef?: PlRef) {
  app.model.data.countsRef = inputRef;
  if (inputRef)
    app.model.data.title = app.model.outputs.countsOptions?.find((o) =>
      plRefsEqual(o.ref, inputRef),
    )?.label;
  else app.model.data.title = undefined;
}

const covariateOptions = computed(() => {
  return (
    app.model.outputs.metadataOptions?.map((v) => ({
      value: v.ref,
      label: v.label,
    })) ?? []
  );
});

const contrastFactorOptions = computed(() => {
  return app.model.data.covariateRefs.map((ref) => ({
    value: ref,
    label: covariateOptions.value.find((m) => m.value.name === ref.name)?.label ?? "",
  }));
});

// Get all possible numerator/denominator values
const uniqueOptions = useWatchFetch(
  () => app.model.outputs.numeratorOptions,
  async (pframeHandle) => {
    if (!pframeHandle) {
      return undefined;
    }
    // Get ID of first pcolumn in the pframe (the only one we will access)
    const pFrame = new PFrameImpl(pframeHandle);
    const list = await pFrame.listColumns();
    const id = list?.[0].columnId;
    if (!id) {
      return undefined;
    }
    // Get unique values of that first pcolumn
    const response = await pFrame.getUniqueValues({ columnId: id, filters: [], limit: 1000000 });
    if (!response) {
      return undefined;
    }
    return [...response.values.data].map((v) => ({ value: String(v), label: String(v) }));
  },
);

// Groups are values of the contrast factor, so a new factor makes the old
// choices meaningless. The reset happens here, on the user's own change, and
// not in a watcher on the fetched values.
function setContrastFactor(ref?: PlRef) {
  app.model.data.contrastFactor = ref;
  app.model.data.numerators = [];
  app.model.data.denominator = undefined;
}

// Chosen groups that the current contrast factor no longer has, for example
// after the upstream data changed.
const missingGroups = computed(() => {
  const available = uniqueOptions.value?.map((op) => op.value);
  if (available === undefined) return [];
  const chosen = [...app.model.data.numerators];
  if (app.model.data.denominator !== undefined) chosen.push(app.model.data.denominator);
  return chosen.filter((group) => !available.includes(group));
});

// Get all possible denominator values (excluding numerators)
const denominatorOptions = computed(() => {
  return uniqueOptions.value?.filter((op) => !app.model.data.numerators.includes(op.value));
});

// Get error logs
const errorLogs = useWatchFetch(
  () => app.model.outputs.deErrors,
  async (pframeHandle) => {
    if (!pframeHandle) {
      return undefined;
    }
    // Get ID of first pcolumn in the pframe (the only one we will access)
    const pFrame = new PFrameImpl(pframeHandle);
    const list = await pFrame.listColumns();
    const id = list?.[0].columnId;
    if (!id) {
      return undefined;
    }
    // Get unique values of that first pcolumn
    const response = await pFrame.getUniqueValues({ columnId: id, filters: [], limit: 1000000 });
    if (!response) {
      return undefined;
    }
    if (response.values.data.length === 0) {
      return undefined;
    }
    return response.values.data.join("\n");
  },
);
</script>

<template>
  <PlBlockPage>
    <template #title>Single-Cell Differential Expression Analysis</template>
    <template #append>
      <PlBtnGhost @click.stop="() => (data.settingsOpen = true)">
        Settings
        <template #append>
          <PlMaskIcon24 name="settings" />
        </template>
      </PlBtnGhost>
    </template>
    <PlAlert v-if="errorLogs.value !== undefined" type="warn" icon>
      {{ errorLogs.value }}
    </PlAlert>
    <PlAgDataTableV2
      v-model="app.model.data.tableState"
      :settings="tableSettings"
      show-columns-panel
      show-export-button
      no-rows-text="There are no differentially expressed genes between the chosen conditions"
    />
    <PlSlideModal v-model="data.settingsOpen">
      <template #title>Settings</template>
      <PlDropdownRef
        v-model="app.model.data.countsRef"
        :options="app.model.outputs.countsOptions"
        label="Select dataset"
        required
        clearable
        @update:model-value="setInput"
      />
      <PlDropdownMulti
        v-model="app.model.data.covariateRefs"
        :options="covariateOptions"
        label="Covariates"
        required
      >
        <template #tooltip>
          <div>
            <strong>Covariates for differential expression</strong><br />
            Select metadata variables to include in the differential expression model. Covariates
            help control for confounding factors (e.g., batch effects, cell cycle, technical
            variables) and improve the accuracy of differential expression results.
          </div>
        </template>
      </PlDropdownMulti>
      <PlDropdown
        :model-value="app.model.data.contrastFactor"
        :options="contrastFactorOptions"
        label="Contrast factor"
        required
        @update:model-value="setContrastFactor"
      >
        <template #tooltip>
          <div>
            <strong>Primary factor for comparison</strong><br />
            Choose the main categorical variable for differential expression analysis. This defines
            the groups you want to compare (e.g., cell types, treatment conditions, disease states).
          </div>
        </template>
      </PlDropdown>
      <PlDropdownMulti
        v-model="app.model.data.numerators"
        :options="uniqueOptions.value"
        label="Numerator (test group)"
        required
      >
        <template #tooltip>
          <div>
            <strong>Test condition or group</strong><br />
            Select the experimental condition or group of interest. Genes upregulated in this group
            will have positive log2 fold-changes in the results.
          </div>
        </template>
      </PlDropdownMulti>
      <PlDropdown
        v-model="app.model.data.denominator"
        :options="denominatorOptions"
        label="Denominator (reference group)"
        required
      >
        <template #tooltip>
          <div>
            <strong>Reference or control group</strong><br />
            Select the baseline condition or control group for comparison. Expression levels in the
            numerator group will be compared relative to this reference group.
          </div>
        </template>
      </PlDropdown>

      <PlAlert v-if="missingGroups.length > 0" type="warn">
        {{ `These groups are not present in the contrast factor: ${missingGroups.join(", ")}` }}
      </PlAlert>

      <PlAccordionSection label="SIGNIFICANCE THRESHOLDS">
        <PlRow>
          <PlNumberField
            v-model="app.model.data.log2FCThreshold"
            label="Log2(FC)"
            :minValue="0"
            :step="0.1"
          >
            <template #tooltip>
              <div>
                <strong>Log2 fold-change threshold</strong><br />
                Minimum absolute log2 fold-change required for a gene to be considered
                differentially expressed. Higher values (e.g., 1.0-2.0) identify more dramatically
                changed genes, while lower values (e.g., 0.5-0.8) capture subtler expression
                differences.
              </div>
            </template>
          </PlNumberField>
          <PlNumberField
            v-model="app.model.data.pAdjThreshold"
            label="Adjusted p-value"
            :minValue="0"
            :maxValue="1"
            :step="0.01"
          >
            <template #tooltip>
              <div>
                <strong>Statistical significance threshold</strong><br />
                Maximum adjusted p-value for differential expression significance. Standard values
                are 0.05 (5% false discovery rate) or 0.01 (1% false discovery rate) for more
                stringent filtering.
              </div>
            </template>
          </PlNumberField>
        </PlRow>
        <!-- Add warnings if selected threshold are out of most commonly used bounds -->
        <PlAlert v-if="app.model.data.pAdjThreshold > 0.05" type="warn">
          {{
            "Warning: The selected adjusted p-value threshold is higher than the commonly recommended 0.05"
          }}
        </PlAlert>
        <PlAlert v-if="app.model.data.log2FCThreshold < 0.6" type="warn">
          {{
            "Warning: The selected Log2(FC) threshold may be too low for identifying robust differentially expressed genes"
          }}
        </PlAlert>
      </PlAccordionSection>
    </PlSlideModal>
  </PlBlockPage>
</template>
