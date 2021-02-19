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
    pagination: { pageIndex: 0, pageSize: 10 },
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
    pagination: { pageIndex: 0, pageSize: 10 },
  },
  rp: {
    search: '',
    hidden: ['First author', 'Last Author', 'Pubmed ID/BioRxiv', 'DOI'],
    pagination: { pageIndex: 0, pageSize: 10 },
  },
  cm: { search: '', hidden: ['Programming'] },
  pagination: { pageIndex: 0, pageSize: 10 },
});

export const { actions, reducer } = createSlice({
  name: 'publications',
  initialState: getInitialState(),
  reducers: {
    mergeState: mergeObject,
    resetModal: getInitialState,
  },
});
