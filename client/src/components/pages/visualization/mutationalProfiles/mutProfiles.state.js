import { atom, selector } from "recoil";
import axios from "axios";
import { visualizationState } from "../visualization.state";
import SBS96 from "../../../controls/plotly/mutationalSignature/sbs96";
import DBS78 from "../../../controls/plotly/mutationalSignature/dbs78";
import ID83 from "../../../controls/plotly/mutationalSignature/id83";

export const defaultFormState = {
  profile: {},
  // filtered: [],
  // selectName: '',
  // selectProfile: '',
  // selectMatrix: '',
  // selectFilter: '',
};

export const formState = atom({
  key: "visualization.mutProfiles.formState",
  default: defaultFormState,
});

// const defaultOptions = {
//   nameOptions: [],
//   profileOptions: [],
//   matrixOptions: [],
//   filterOptions: [],
// };

// export const formOptions = selector({
//   key: 'visualization.mutProfiles.formOptions',
//   get: async ({ get }) => {
//     const { svgList } = get(visualizationState);
//     try {
//       const { data } = await axios.post('web/visualizationWrapper', {
//         fn: 'getTreeLeaf',
//       });

//       return data.output;
//       return [];
//     } catch (error) {
//       return null;
//     }
//   },
// });

export const getPlot = selector({
  key: "visualization.mutProfiles.plot",
  get: async ({ get }) => {
    try {
      const { profile } = get(formState);

      if (Object.keys(profile).length) {
        const { data } = await axios.post("web/querySignature", {
          filter: { ...profile.value }, // filter given an object of key:values you want
          properties: ["MutationType", "Contribution"], // return objects containing these properties
        });

        console.log(profile.value.Profile);

        if (profile.value.Profile === "SBS96") {
          console.log(data);
          const { traces, layout } = SBS96(data);
          return { data: [...traces], layout, config: {} };
        } else if (profile.value.Profile === "DBS78") {
          console.log(data);
          const { traces, layout } = DBS78(data);
          return { data: [...traces], layout, config: {} };
        } else if (profile.value.Profile === "ID83") {
          console.log(data);
          const { traces, layout } = ID83(data);
          return { data: [...traces], layout, config: {} };
        } else return { data: [], layout: {}, config: {} };
      } else return { data: [], layout: {}, config: {} };
    } catch (error) {
      return { data: [], layout: {}, config: {} };
    }
  },
});
