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

/** Leaf color options available for user-uploaded data (seqmatrix-only). */
export const userColorOptions = [
  {
    label: 'Dominant Mutation',
    value: 'Dmut',
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
  signatureSetName: 'COSMIC_v3_Signatures_GRCh37_SBS96',
  profile: 'SBS',
  matrix: 96,
};

export const defaultUserFormState = {
  cancerType: null,
  showLabels: false,
  color: userColorOptions[0],
  signatureSetName: null,
  profile: 'SBS',
  matrix: 96,
};

export const formState = atom({
  key: 'treeLeaf.formState',
  default: defaultFormState,
});

export const defaultTreeLeafData = { links: [], nodes: [] };

export const treeLeafDataState = atom({
  key: 'treeLeaf.coordinateState',
  default: defaultTreeLeafData,
});

export const graphDataSelector = selectorFamily({
  key: 'treeLeaf.plotData',
  get:
    (params) =>
    async ({ get }) => {
      try {
        if (params?.source === 'user' && !params?.userId) {
          return null;
        }
        const response = await axios.post('api/treeLeaf', params);
        const data = response.data.output;
        if (data?.error || data?.uncaughtError) {
          return {
            error:
              data.error ||
              data.uncaughtError ||
              'Tree and Leaf calculation failed.',
          };
        }
        return data;
      } catch (error) {
        console.error(error);
        const message =
          error?.response?.data?.error ||
          error?.response?.data?.message ||
          error?.message ||
          'Failed to load Tree and Leaf data.';
        return { error: message };
      }
    },
});
