import { atom, selector } from 'recoil';
import axios from 'axios';
import sigPatternData from './sigPatternData.json';

export const colorOptions = [
  {
    label: 'Cosine Similarity',
    value: 'Cosine_similarity',
    continuous: true,
  },
  {
    label: 'Dominant Mutation',
    value: 'Dmut',
    continuous: false,
  },
  {
    label: 'Dominant Signature',
    value: 'Dsig',
    continuous: false,
  },
  {
    label: 'Cancer Type',
    value: 'Cancer_Type',
    continuous: false,
  },
];

export const defaultFormState = {
  showLabels: false,
  color: {
    label: 'Cosine Similarity',
    value: 'Cosine_similarity',
    continuous: true,
  },
};

export const formState = atom({
  key: 'treeLeaf.formState',
  default: defaultFormState,
});

export const defaultTreeLeafData = { links: [], nodes: [] };

export const graphDataSelector = selector({
  key: 'treeLeaf.plotData',
  get: async ({ get }) => {
    try {
      // const { data } = await axios.post('web/visualizationWrapper', {
      //   fn: 'getTreeLeaf',
      // });

      // return data.output;
      return sigPatternData;
    } catch (error) {
      return null;
    }
  },
});
