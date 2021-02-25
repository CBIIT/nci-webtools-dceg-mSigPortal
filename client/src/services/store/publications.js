import { createSlice } from '@reduxjs/toolkit';
import { mergeObject } from './utils';

// data is loaded in components/app.js
export const getInitialState = () => ({
  orA: {
    globalFilter: '',
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
    globalFilter: '',
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
    globalFilter: '',
    hidden: ['First author', 'Last Author', 'Pubmed ID/BioRxiv', 'DOI'],
    pagination: { pageIndex: 0, pageSize: 10 },
  },
  cm: {
    globalFilter: '',
    hidden: [
      'Programming',
      'First author',
      'Last Author',
      'Pubmed ID/BioRxiv',
      'DOI',
    ],
    pagination: { pageIndex: 0, pageSize: 10 },
  },
});

export const { actions, reducer } = createSlice({
  name: 'publications',
  initialState: getInitialState(),
  reducers: {
    mergeState: mergeObject,
    resetModal: getInitialState,
  },
});
