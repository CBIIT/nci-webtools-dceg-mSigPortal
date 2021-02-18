import { createSlice } from '@reduxjs/toolkit';
import { mergeObject } from './utils';

// data is loaded in components/app.js
export const getInitialState = () => ({
  orA: {
    hidden: [
      'Disease/Phenotype/Exposure',
      'First author',
      'Last Author',
      'Pubmed ID/BioRxiv',
      'DOI',
    ],
  },
  orB: {
    hidden: [
      'Disease/Phenotype/Exposure',
      'First author',
      'Last Author',
      'Pubmed ID/BioRxiv',
      'DOI',
    ],
  },
  rp: { hidden: ['First author', 'Last Author', 'Pubmed ID/BioRxiv', 'DOI'] },
  cm: { hidden: ['Programming'] },
});

export const { actions, reducer } = createSlice({
  name: 'publications',
  initialState: getInitialState(),
  reducers: {
    mergeState: mergeObject,
    resetModal: getInitialState,
  },
});
