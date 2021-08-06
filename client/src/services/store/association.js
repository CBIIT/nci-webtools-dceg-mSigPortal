import { createSlice } from '@reduxjs/toolkit';
import { mergeObject } from './utils';

export const getInitialState = () => ({
  associationState: {
    displayTab: 'univariate',
    openSidebar: true,
    loadingData: false,
    submitted: false,
    exposureSignature: [],
    assocVarData: [],
    expVarList: [],
    source: 'public',
    study: '',
    studyOptions: [],
    strategy: '',
    strategyOptions: [],
    rsSet: '',
    rsSetOptions: [],
    cancer: '',
    cancerOptions: [],
    assocTable: {
      hidden: [],
      pagination: {
        pageIndex: 0,
        pageSize: 5,
      },
    },
  },
  univariate: {
    loadingParams: false,
    loadingCalculate: false,
    loadingRecalculate: false,
    submitted: false,
    error: false,
    loadedParameters: false,
    assocVariant: '',
    expVariant: '',
    projectID: '',
    plotPath: '',
    dataPath: '',
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
    resultsTable: {
      data: [],
      hidden: [],
      pagination: {
        pageIndex: 0,
        pageSize: 5,
      },
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
