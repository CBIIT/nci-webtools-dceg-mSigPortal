import { createSlice } from '@reduxjs/toolkit';
import { mergeObject } from './utils';

// data is loaded in components/app.js
export const getInitialState = () => ({
  orA: {
    search: '',
    hidden: [
      'Disease/Phenotype/Exposure',
      'First author',
      'Last Author',
      'Pubmed ID/BioRxiv',
      'DOI',
    ],
  },
  orB: {
    search: '',
    hidden: [
      'Disease/Phenotype/Exposure',
      'First author',
      'Last Author',
      'Pubmed ID/BioRxiv',
      'DOI',
    ],
  },
  rp: {
    search: '',
    hidden: ['First author', 'Last Author', 'Pubmed ID/BioRxiv', 'DOI'],
  },
  cm: { search: '', hidden: ['Programming'] },
});

export const { actions, reducer } = createSlice({
  name: 'publications',
  initialState: getInitialState(),
  reducers: {
    mergeState: mergeObject,
    resetModal: getInitialState,
  },
});
