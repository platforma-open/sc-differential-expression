import { model } from '@platforma-open/milaboratories.sc-differential-expression.model';
import { defineApp } from '@platforma-sdk/ui-vue';
import MainPage from './pages/MainPage.vue';
import VolcanoPage from './pages/GraphPage.vue';

export const sdkPlugin = defineApp(model, () => {
  return {
    routes: {
      '/': () => MainPage,
      '/volcano': () => VolcanoPage,
    },
  };
});

export const useApp = sdkPlugin.useApp;
