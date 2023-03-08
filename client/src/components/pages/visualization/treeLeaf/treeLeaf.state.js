import { atom, selectorFamily } from 'recoil';
import axios from 'axios';

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
  cancerType: null,
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

export const graphDataSelector = selectorFamily({
  key: 'treeLeaf.plotData',
  get: (params) => async ({ get }) => {
    try {
      const response = await axios.post('web/treeLeaf', params);
      const data = response.data.output;
      console.log(data);
      return data;
    } catch (error) {
      console.error(error);
      return null;
    }
  },
});
