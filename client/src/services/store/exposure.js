import { createSlice } from '@reduxjs/toolkit';
import { mergeObject } from './utils';

export const getInitialState = () => ({
  publicForm: {
    study: '',
    strategy: '',
    signatureSet: '',
    cancer: '',
    useAllCancer: false,
  },
  userForm: {},
  main: {
    displayTab: 'instructions',
    source: 'public',
    loading: false,
    projectID: '',
    openSidebar: true,
    submitted: false,
  },

  msAssociation: {
    both: true,
    signatureName1: '',
    signatureName2: '',
    plotPath: '',
    txtPath: '',
    debugR: '',
    err: '',
    loading: false,
  },
  msLandscape: {
    variableFile: '',
    plotPath: '',
    txtPath: '',
    debugR: '',
    err: '',
    loading: false,
  },
  msIndividual: {
    sample: '',
    plotPath: '',
    debugR: '',
    err: '',
    loading: false,
  },
});

export const { actions, reducer } = createSlice({
  name: 'exposure',
  initialState: getInitialState(),
  reducers: {
    mergeExposure: mergeObject,
    resetExposure: getInitialState,
  },
});
