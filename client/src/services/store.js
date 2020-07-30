import { createSlice, configureStore, combineReducers } from '@reduxjs/toolkit';

export const getInitialState = () => ({
  visualize: {
    inputFormat: 'vcf',
    selectedGenome: 'GRCh37',
    experimentalStrategy: 'WGS',
    mutationSplit: 'False',
    collapseSample: 'False',
    mutationFilter: '',
    queueMode: false,
    email: '',
    openSidebar: true,
    disableParameters: false,
    storeFile: '',
    submitted: false,
    exampleData: 'assets/exampleInput/demo_input_multi.vcf.gz',
    loading: {
      active: false,
      content: null,
      showIndicator: false,
    },
  },
  visualizeResults: {
    error: '',
    projectID: '',
    displayTab: 'python',
    pyPlotURL: '',
    csPlotURL: '',
    pyTab: {
      mapping: [],
      filtered: [],
      selectName: '',
      selectProfile: '',
      selectMatrix: '',
      selectTag: '',
      nameOptions: [],
      profileOptions: [],
      matrixOptions: [],
      tagOptions: [],
      debugPy: { stdout: '', stderr: '' },
    },
    cosineSimilarity: {
      profileType1: '',
      matrixSize: '',
      matrixOptions: [],
      profileType2: '',
      signatureSet: '',
      signatureSetOptions: [],
      // selectProfile: '',
      // selectSigFormula: 'signature',
      // sigFormula: 'SBS5',
      rPlots: [],
      rPlotIndex: '',
      results: [],
      debugR: [],
      submitOverlay: false,
      refSigOverlay: false,
    },
  },
  error: {
    visible: false,
    message: `An error occured when requesting data. If this problem persists, please contact the administrator at <a href="mailto:mSigPortalWebAdmin@cancer.gov">mSigPortalWebAdmin@cancer.gov</a>.`,
  },
});

const visualizeSlice = createSlice({
  name: 'visualize',
  initialState: getInitialState().visualize,
  reducers: {
    updateVisualize: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});

const visualizeResultsSlice = createSlice({
  name: 'visualizeResults',
  initialState: getInitialState().visualizeResults,
  reducers: {
    updateVisualizeResults: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});

const errorSlice = createSlice({
  name: 'error',
  initialState: getInitialState().error,
  reducers: {
    updateError: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});

const rootReducer = combineReducers({
  visualize: visualizeSlice.reducer,
  visualizeResults: visualizeResultsSlice.reducer,
  error: errorSlice.reducer,
});

export const store = configureStore({
  reducer: rootReducer,
  preloadedState: getInitialState(),
});

export const { updateVisualize } = visualizeSlice.actions;

export const { updateVisualizeResults } = visualizeResultsSlice.actions;

export const { updateError } = errorSlice.actions;

export function dispatchVisualize(obj) {
  store.dispatch(updateVisualize(obj));
}

export function dispatchVisualizeResults(obj) {
  store.dispatch(updateVisualizeResults(obj));
}

export function dispatchError(msg) {
  store.dispatch(
    updateError({
      visible: true,
      message: msg,
    })
  );
}
