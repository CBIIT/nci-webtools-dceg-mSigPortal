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

export const SIGNATURE_SET_NAME = 'COSMIC_v3_Signatures_GRCh37_SBS96';
export const PROFILE = 'SBS';
export const MATRIX = 96;

/**
 * Builds the initial Tree and Leaf form state synchronously from the sidebar form
 */
export function getInitialFormState({ isUser, publicForm }) {
  if (isUser) {
    return { color: userColorOptions[0], searchSamples: [] };
  }
  const value = publicForm?.cancer?.value;
  const known = (publicForm?.cancers || []).some(
    (c) => c.value === value && c.value !== '*ALL'
  );
  return {
    color: colorOptions[0],
    searchSamples: [],
    cancer: known ? value : '',
  };
}

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
        if (data?.error) {
          return {
            error: data.error,
          };
        } else if (data?.uncaughtError) {
          return {
            error: `An error occurred with the selected study: ${data.uncaughtError}`,
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
