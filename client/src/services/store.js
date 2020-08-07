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
    downloads: [],
    summary: [],
    matrixList: [],
    statistics: '',
  },
  pyTab: {
    filtered: [],
    selectName: '',
    selectProfile: '',
    selectMatrix: '',
    selectTag: '',
    nameOptions: [],
    profileOptions: [],
    matrixOptions: [],
    tagOptions: [],
    plotURL: '',
    debug: { stdout: '', stderr: '' },
  },
  cosineSimilarity: {
    withinProfileType: '',
    withinMatrixSize: '',
    withinMatrixOptions: [],
    refProfileType2: '',
    refSignatureSet: '',
    refSignatureSetOptions: [],
    withinPlotPath: '',
    withinTxtPath: '',
    refPlotPath: '',
    refTxtPath: '',
    withinPlotURL: '',
    refPlotURL: '',
    displayWithin: true,
    displayRefSig: true,
    debugR: [],
    withinSubmitOverlay: false,
    refSubmitOverlay: false,
  },
  profileComparison: {
    withinProfileType: '',
    withinSampleName1: '',
    withinSampleName2: '',
    refProfileType: '',
    refSampleName: '',
    refSignatureSet: '',
    refSignatureSetOptions: [],
    refCompare: 'SBS5',
    withinPlotPath: '',
    refPlotPath: '',
    withinPlotURL: '',
    refPlotURL: '',
    displayWithin: true,
    displayRefSig: true,
    debugR: [],
    withinSubmitOverlay: false,
    refSubmitOverlay: false,
  },
  pca: {
    profileType: '',
    signatureSet: '',
    signatureSetOptions: [],
    eig: '',
    pca1: '',
    pca2: '',
    heatmap: '',
    pca1Data: '',
    pca2Data: '',
    heatmapData: '',
    eigURL: '',
    pca1URL: '',
    pca2URL: '',
    heatmapURL: '',
    displayPCA: true,
    debugR: [],
    submitOverlay: false,
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

const pyTabSlice = createSlice({
  name: 'pyTab',
  initialState: getInitialState().pyTab,
  reducers: {
    updatePyTab: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});

const cosineSimilaritySlice = createSlice({
  name: 'cosineSimilarity',
  initialState: getInitialState().cosineSimilarity,
  reducers: {
    updateCosineSimilarity: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});

const profileComparisonSlice = createSlice({
  name: 'profileComparison',
  initialState: getInitialState().profileComparison,
  reducers: {
    updateProfileComparison: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});

const pcaSlice = createSlice({
  name: 'pca',
  initialState: getInitialState().pca,
  reducers: {
    updatePCA: (state, action) => {
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
  pyTab: pyTabSlice.reducer,
  cosineSimilarity: cosineSimilaritySlice.reducer,
  profileComparison: profileComparisonSlice.reducer,
  pca: pcaSlice.reducer,
  error: errorSlice.reducer,
});

export const store = configureStore({
  reducer: rootReducer,
  preloadedState: getInitialState(),
});

export const { updateVisualize } = visualizeSlice.actions;
export const { updateVisualizeResults } = visualizeResultsSlice.actions;
export const { updatePyTab } = pyTabSlice.actions;
export const { updateCosineSimilarity } = cosineSimilaritySlice.actions;
export const { updateProfileComparison } = profileComparisonSlice.actions;
export const { updatePCA } = pcaSlice.actions;
export const { updateError } = errorSlice.actions;

export function dispatchVisualize(obj) {
  store.dispatch(updateVisualize(obj));
}

export function dispatchVisualizeResults(obj) {
  store.dispatch(updateVisualizeResults(obj));
}

export function dispatchPyTab(obj) {
  store.dispatch(updatePyTab(obj));
}

export function dispatchCosineSimilarity(obj) {
  store.dispatch(updateCosineSimilarity(obj));
}

export function dispatchProfileComparison(obj) {
  store.dispatch(updateProfileComparison(obj));
}

export function dispatchPCA(obj) {
  store.dispatch(updatePCA(obj));
}

export function dispatchError(msg) {
  store.dispatch(
    updateError({
      visible: true,
      message: msg,
    })
  );
}
