import { createSlice } from '@reduxjs/toolkit';
import { mergeObject } from './utils';

export const getInitialState = () => ({
  openSidebar: true,
  loading: false,
  submitted: false,
  err: false,
  loadedParameters: false,
  exposureSignature: [],
  assocVarData: [],
  expVarList: [],
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
  regression: false,
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
  hidden: [],
  pagination: {
    pageIndex: 0,
    pageSize: 10,
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
