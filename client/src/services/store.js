import { createSlice, configureStore, combineReducers } from '@reduxjs/toolkit';

export const getInitialState = () => ({
  visualize: {
    source: 'user',
    inputFormat: 'vcf',
    selectedGenome: 'GRCh37',
    experimentalStrategy: 'WGS',
    mutationSplit: 'False',
    collapseSample: 'False',
    mutationFilter: '',
    queueMode: false,
    email: '',
    study: '',
    studyOptions: [],
    cancerType: '',
    cancerTypeOptions: [],
    pubExperimentalStrategy: '',
    pubExperimentOptions: [],
    pDataOptions: [],
    openSidebar: true,
    storeFilename: '',
    bedFilename: '',
    submitted: false,
    exampleData: 'assets/exampleInput/demo_input_multi.vcf.gz',
    bedData: 'assets/exampleInput/demo_input_bed.bed',
    loadingPublic: false,
    loading: {
      active: false,
      content: null,
      showIndicator: false,
    },
  },
  visualizeResults: {
    error: '',
    projectID: '',
    // projectID: '04640afa-5530-4790-91db-6ff695a51478',
    displayTab: 'profilerSummary',
    downloads: [],
    svgList: [],
    matrixList: [],
    statistics: '',
  },
  profilerSummary: {
    plotPath: '',
    plotURL: '',
    err: '',
    debugR: '',
    displayDebug: false,
    loading: false,
  },
  mutationalProfiles: {
    filtered: [],
    selectName: '',
    selectProfile: '',
    selectMatrix: '',
    selectFilter: '',
    nameOptions: [],
    profileOptions: [],
    matrixOptions: [],
    filterOptions: [],
    plotURL: '',
    debug: { stdout: '', stderr: false },
    displayDebug: false,
  },
  cosineSimilarity: {
    withinProfileType: '',
    withinMatrixSize: '',
    withinMatrixOptions: [],
    refProfileType: '',
    refSignatureSet: '',
    refSignatureSetOptions: [],
    userProfileType: '',
    userMatrixSize: '',
    userMatrixOptions: [],
    pubStudy: '',
    pubCancerType: '',
    pubCancerTypeOptions: [],
    withinPlotPath: '',
    withinTxtPath: '',
    refPlotPath: '',
    refTxtPath: '',
    pubPlotPath: '',
    pubTxtPath: '',
    withinPlotURL: '',
    refPlotURL: '',
    pubPlotURL: '',
    displayWithin: true,
    displayRefSig: true,
    displayPublic: true,
    withinErr: false,
    refErr: false,
    pubErr: false,
    debugR: [],
    displayDebug: false,
    withinSubmitOverlay: false,
    refSubmitOverlay: false,
    pubSubmitOverlay: false,
  },
  mutationalPattern: {
    proportion: '',
    pattern: '',
    txtPath: '',
    plotPath: '',
    plotURL: '',
    barPath: '',
    barURL: '',
    display: true,
    err: false,
    debugR: [],
    displayDebug: false,
    submitOverlay: false,
  },
  profileComparison: {
    withinProfileType: '',
    withinSampleName1: '',
    withinSampleName2: '',
    refProfileType: '',
    refSampleName: '',
    refSignatureSet: '',
    refSignatureSetOptions: [],
    refSignatures: [],
    refCompare: '',
    userProfileType: '',
    userMatrixSize: '',
    userMatrixOptions: '',
    userSampleName: '',
    pubSampleName: '',
    pubSampleOptions: [],
    pubStudy: '',
    pubCancerType: '',
    pubCancerTypeOptions: [],
    withinPlotPath: '',
    refPlotPath: '',
    withinPlotURL: '',
    refPlotURL: '',
    pubPlotPath: '',
    pubPlotURL: '',
    displayWithin: true,
    displayRefSig: true,
    withinErr: false,
    refErr: false,
    pubErr: false,
    debugR: [],
    displayDebug: false,
    withinSubmitOverlay: false,
    refSubmitOverlay: false,
    pubSubmitOverlay: false,
  },
  pca: {
    profileType: '',
    signatureSet: '',
    signatureSetOptions: [],
    pca1: '',
    pca2: '',
    pca3: '',
    heatmap: '',
    pca2Data: '',
    pca3Data: '',
    heatmapData: '',
    pca1URL: '',
    pca2URL: '',
    pca3URL: '',
    heatmapURL: '',
    displayPCA: true,
    pca1Err: false,
    pca2Err: false,
    pca3Err: false,
    heatmapErr: false,
    displayPCA: true,
    debugR: [],
    displayDebug: false,
    submitOverlay: false,
    userProfileType: '',
    userMatrixSize: '',
    userMatrixOptions: [],
    pubStudy: '',
    pubCancerType: '',
    pubCancerTypeOptions: [],
    pubPca1: '',
    pubPca2: '',
    pubPca3: '',
    pubPca2Data: '',
    pubPca3Data: '',
    pubPca1URL: '',
    pubPca2URL: '',
    pubPca3URL: '',
    displayPub: true,
    pubPca1Err: false,
    pubPca2Err: false,
    pubPca3Err: false,
    pubSubmitOverlay: false,
  },

  exploring: {
    displayTab: 'signatureExploring',
    signatureAccordion: {
      referenceSignatures: true,
      mutationalSignatureProfile: true,
      cosineSimilarity: true,
      mutationalSignatureComparison: true,
    },
    exposureAccordion: {},
    refSigData: {},
    projectID: '',
    source: 'public',
    openSidebar: false,
    submitted: false,
  },
  expRefSig: {
    plotPath: '',
    plotURL: 'assets/images/mSigPortalReferenceSignatures.svg',
    debugR: '',
    err: '',
    displayDebug: false,
    loading: false,
  },
  expMutationalProfiles: {
    signatureSource: '',
    signatureSourceOptions: [],
    profileName: '',
    profileNameOptions: [],
    refSignatureSet: '',
    refSignatureSetOptions: [],
    strategy: '',
    strategyOptions: [],
    signatureName: '',
    signatureNameOptions: [],

    plotPath: '',
    plotURL: '',
    debugR: '',
    err: '',
    displayDebug: false,
    loading: false,
  },
  expCosineSimilarity: {
    profileName: '',
    profileNameOptions: [],
    refSignatureSet1: '',
    refSignatureSet2: '',
    refSignatureSetOptions1: [],
    refSignatureSetOptions2: [],
    plotPath: '',
    plotURL: '',
    txtPath: '',
    debugR: '',
    err: '',
    displayDebug: false,
    loading: false,
  },
  expMutationalSigComparison: {
    profileName: '',
    profileNameOptions: [],
    refSignatureSet1: '',
    refSignatureSetOptions1: [],
    refSignatureSet2: '',
    refSignatureSetOptions2: [],
    signatureName1: '',
    signatureNameOptions1: [],
    signatureName2: '',
    signatureNameOptions2: [],
    plotPath: '',
    plotURL: '',
    debugR: '',
    err: '',
    displayDebug: false,
    loading: false,
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

const profilerSummarySlice = createSlice({
  name: 'profilerSummary',
  initialState: getInitialState().profilerSummary,
  reducers: {
    updateProfilerSummary: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});

const mutationalProfilesSlice = createSlice({
  name: 'mutationalProfiles',
  initialState: getInitialState().mutationalProfiles,
  reducers: {
    updateMutationalProfiles: (state, action) => {
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

const mutationalPatternSlice = createSlice({
  name: 'mutationalPattern',
  initialState: getInitialState().mutationalPattern,
  reducers: {
    updateMutationalPattern: (state, action) => {
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

const exploringSlice = createSlice({
  name: 'exploring',
  initialState: getInitialState().exploring,
  reducers: {
    updateExploring: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});

const expRefSigSlice = createSlice({
  name: 'expRefSig',
  initialState: getInitialState().expRefSig,
  reducers: {
    updateExpRefSig: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});

const expMutationalProfilesSlice = createSlice({
  name: 'expMutationalProfiles',
  initialState: getInitialState().expMutationalProfiles,
  reducers: {
    updateExpMutationalProfiles: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});

const expCosineSimilaritySlice = createSlice({
  name: 'expCosineSimilarity',
  initialState: getInitialState().expCosineSimilarity,
  reducers: {
    updateExpCosineSimilarity: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});

const expMutationalSigComparisonSlice = createSlice({
  name: 'expMutationalSigComparison',
  initialState: getInitialState().expMutationalSigComparison,
  reducers: {
    updateExpMutationalSigComparison: (state, action) => {
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
  profilerSummary: profilerSummarySlice.reducer,
  mutationalProfiles: mutationalProfilesSlice.reducer,
  cosineSimilarity: cosineSimilaritySlice.reducer,
  mutationalPattern: mutationalPatternSlice.reducer,
  profileComparison: profileComparisonSlice.reducer,
  pca: pcaSlice.reducer,
  error: errorSlice.reducer,
  exploring: exploringSlice.reducer,
  expRefSig: expRefSigSlice.reducer,
  expMutationalProfiles: expMutationalProfilesSlice.reducer,
  expCosineSimilarity: expCosineSimilaritySlice.reducer,
  expMutationalSigComparison: expMutationalSigComparisonSlice.reducer,
});

export const store = configureStore({
  reducer: rootReducer,
  preloadedState: getInitialState(),
});

export const { updateVisualize } = visualizeSlice.actions;
export const { updateVisualizeResults } = visualizeResultsSlice.actions;
export const { updateProfilerSummary } = profilerSummarySlice.actions;
export const { updateMutationalProfiles } = mutationalProfilesSlice.actions;
export const { updateCosineSimilarity } = cosineSimilaritySlice.actions;
export const { updateMutationalPattern } = mutationalPatternSlice.actions;
export const { updateProfileComparison } = profileComparisonSlice.actions;
export const { updatePCA } = pcaSlice.actions;
export const { updateError } = errorSlice.actions;
export const { updateExploring } = exploringSlice.actions;
export const { updateExpRefSig } = expRefSigSlice.actions;
export const {
  updateExpMutationalProfiles,
} = expMutationalProfilesSlice.actions;
export const { updateExpCosineSimilarity } = expCosineSimilaritySlice.actions;
export const {
  updateExpMutationalSigComparison,
} = expMutationalSigComparisonSlice.actions;

export function dispatchVisualize(obj) {
  store.dispatch(updateVisualize(obj));
}

export function dispatchVisualizeResults(obj) {
  store.dispatch(updateVisualizeResults(obj));
}

export function dispatchProfilerSummary(obj) {
  store.dispatch(updateProfilerSummary(obj));
}

export function dispatchMutationalProfiles(obj) {
  store.dispatch(updateMutationalProfiles(obj));
}

export function dispatchCosineSimilarity(obj) {
  store.dispatch(updateCosineSimilarity(obj));
}

export function dispatchMutationalPattern(obj) {
  store.dispatch(updateMutationalPattern(obj));
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

export function dispatchExploring(obj) {
  store.dispatch(updateExploring(obj));
}

export function dispatchExpRefSig(obj) {
  store.dispatch(updateExpRefSig(obj));
}

export function dispatchExpMutationalProfiles(obj) {
  store.dispatch(updateExpMutationalProfiles(obj));
}

export function dispatchExpCosineSimilarity(obj) {
  store.dispatch(updateExpCosineSimilarity(obj));
}

export function dispatchExpMutationalSigComparison(obj) {
  store.dispatch(updateExpMutationalSigComparison(obj));
}
