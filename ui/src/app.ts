import { platforma } from "@platforma-open/milaboratories.sc-differential-expression.model";
import { defineAppV3 } from "@platforma-sdk/ui-vue";
import MainPage from "./pages/MainPage.vue";
import VolcanoPage from "./pages/GraphPage.vue";

export const sdkPlugin = defineAppV3(platforma, (app) => {
  return {
    progress: () => {
      return app.model.outputs.isRunning;
    },
    routes: {
      "/": () => MainPage,
      "/volcano": () => VolcanoPage,
    },
  };
});

export const useApp = sdkPlugin.useApp;
