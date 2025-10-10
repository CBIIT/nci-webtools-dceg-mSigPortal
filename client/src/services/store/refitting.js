import { createSlice } from '@reduxjs/toolkit';
import { mergeObject } from './utils';

export const getInitialState = () => ({
  main: {
    displayTab: 'instructions',
    openSidebar: true,
    submitted: false,
    loading: false,
    id: null,
  },
  userForm: {
    signatureType: 'SBS',
    referenceGenome: 'hg19',
    mafFile: null,
    genomicFile: null,
    clinicalFile: null,
    email: '',
    jobName: '',
    matchOnOncotree: false,
    outputFilename: '',
  },
});

export const { actions, reducer } = createSlice({
  name: 'refitting',
  initialState: getInitialState(),
  reducers: {
    mergeRefitting: mergeObject,
    resetRefitting: getInitialState,
  },
});