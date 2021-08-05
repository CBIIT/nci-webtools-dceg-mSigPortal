import { createSlice } from '@reduxjs/toolkit';
import { mergeObject } from './utils';

export const getInitialState = () => ({
  openSidebar: true,
  loadingData: false,
  loadingParams: false,
  loadingCalculate: false,
  loadingRecalculate: false,
  submitted: false,
  error: false,
  loadedParameters: false,
  exposureSignature: [],
  assocVarData: [],
  expVarList: [],
  assocVariant: '',
  expVariant: '',
  projectID: '',
  plotPath: '',
  dataPath: '',
  source: 'public',
  study: '',
  studyOptions: [],
  strategy: '',
  strategyOptions: [],
  rsSet: '',
  rsSetOptions: [],
  cancer: '',
  cancerOptions: [],
  dataSource: '',
  dataSourceOptions: [],
  dataType: '',
  dataTypeOptions: [],
  assocVarOptions: [],
  signature: '',
  signatureOptions: [],
  testType: 'nonparametric',
  xlab: '',
  ylab: '',
  variant1: {
    name: '',
    filter: false,
    log2: false,
    collapse: '',
    collapseOptions: [],
  },
  variant2: {
    name: '',
    filter: false,
    log2: false,
    collapse: '',
    collapseOptions: [],
  },
  assocTable: {
    hidden: [],
    pagination: {
      pageIndex: 0,
      pageSize: 5,
    },
  },
  resultsTable: {
    data: [],
    hidden: [],
    pagination: {
      pageIndex: 0,
      pageSize: 5,
    },
  },
});

export const { actions, reducer } = createSlice({
  name: 'association',
  initialState: getInitialState(),
  reducers: {
    mergeAssociation: mergeObject,
    resetAssociation: getInitialState,
  },
});
