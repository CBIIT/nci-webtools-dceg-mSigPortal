import { atom, selector } from 'recoil';
import axios from 'axios';
// import { exposureState } from "../exposure.state";
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

        const { traces, layout } = TMB(data, 'PCAWG');

        return { data: traces, layout: layout, config: {} };
      } else return { data: [], layout: {}, config: {} };
    } catch (error) {
      return { data: [], layout: {}, config: {} };
    }
  },
});
