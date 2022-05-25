import { atom, selector } from 'recoil';
import axios from 'axios';
import sigPatternData from './sigPatternData.json';

export const defaultFormState = {
  showLabels: false,
  color: { label: '', value: '', continuous: false },
};

export const formState = atom({
  key: 'treeLeaf.formState',
  default: defaultFormState,
});

export const defaultTreeLeafData = { links: [], nodes: [] };

export const getGraphData = selector({
  key: 'treeLeaf.plotData',
  get: async ({ get }) => {
    try {
      // const { data } = await axios.post('api/visualizationWrapper', {
      //   fn: 'getTreeLeaf',
      // });

      // return data.output;
      return sigPatternData;
    } catch (error) {
      return null;
    }
  },
});
