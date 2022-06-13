import { atom, selector } from 'recoil';
import axios from 'axios';
import SBS96 from '../../../controls/plotly/mutationalSignature/sbs96';
import DBS78 from '../../../controls/plotly/mutationalSignature/dbs78';
import ID83 from '../../../controls/plotly/mutationalSignature/id83';

export const defaultFormState = {
  option: {},
  // filtered: [],
  // selectName: '',
  // selectProfile: '',
  // selectMatrix: '',
  // selectFilter: '',
};

export const formState = atom({
  key: 'visualization.mutProfiles.formState',
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
  key: 'visualization.mutProfiles.plot',
  get: async ({ get }) => {
    try {
      const { option } = get(formState);

      if (Object.keys(option).length) {
        const { data } = await axios.get('web/querySignature', {
          params: option.value,
        });

        console.log(data);
        if (option.value.profile === 'SBS96') {
          const { traces, layout } = SBS96(data);
          return { data: [...traces], layout, config: {} };
        } else if (option.value.profile === 'DBS78') {
          const { traces, layout } = DBS78(data);
          return { data: [...traces], layout, config: {} };
        } else if (option.value.profile === 'ID83') {
          const { traces, layout } = ID83(data);
          return { data: [...traces], layout, config: {} };
        } else return { data: [], layout: {}, config: {} };
      } else return { data: [], layout: {}, config: {} };
    } catch (error) {
      return { data: [], layout: {}, config: {} };
    }
  },
});
