import { atom, selector } from 'recoil';
import axios from 'axios';
// import { main } from "../exposure.state";
import TMB from '../../../controls/plotly/tmb/tmb';

export const defaultFormState = {
  option: {},
};

export const formState = atom({
  key: 'exposure.mutProfiles.formState',
  default: defaultFormState,
});

export const getPlot = selector({
  key: 'exposure.mutProfiles.plot',
  get: async ({ get }) => {
    try {
      const { option } = get(formState);

      if (Object.keys(option).length) {
        const { data } = await axios.get('web/exposure', {
          params: option.value,
        });

        console.log(option.value.profile);
        if (option.value.profile === 'single') {
        } else {
        }

        const { traces, layout, config } = TMB(data, 'PCAWG');

        return { data: traces, layout: layout, config: config };
      } else return { data: [], layout: {}, config: {} };
    } catch (error) {
      return { data: [], layout: {}, config: {} };
    }
  },
});
