import { createSlice } from '@reduxjs/toolkit';
import { mergeObject } from './utils';

export const getInitialState = () => ({
  inputForm: {
    source: 'public',
    study: '',
    strategy: '',
    cancer: '',
    dataType: '',
    inputFilename: '',
    genome: '',
    exome: '',
  },
  main: {
    displayTab: 'instructions',
    projectID: '',
    openSidebar: true,
    submitted: false,
  },
});

export const { actions, reducer } = createSlice({
  name: 'extraction',
  initialState: getInitialState(),
  reducers: {
    mergeExtraction: mergeObject,
    resetExtraction: getInitialState,
  },
});
