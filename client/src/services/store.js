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
    csWithinURL: '',
    csRefSigURL: '',
    pcWithinURL: '',
    pcRefSigURL: '',
  },
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
    csWithinPlot: '',
    csWithinTxt: '',
    csRefSigPlot: '',
    csRefSigTxt: '',
    displayWithin: true,
    displayRefSig: true,
    debugR: [],
    submitOverlay: false,
  },
  profileComparison: {
    within: {
      profileType: '',
      sampleName1: '',
      sampleName2: '',
    },
    refSig: {
      profileType: '',
      sampleName: '',
      signatureSet: '',
      signatureSetOptions: [],
      compare: 'SBS5',
    },
    withinPlotPath: '',
    refSigPlotPath: '',
    displayWithin: true,
    displayRefSig: true,
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

export function dispatchError(msg) {
  store.dispatch(
    updateError({
      visible: true,
      message: msg,
    })
  );
}
