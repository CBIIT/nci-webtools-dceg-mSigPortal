import { createSlice, configureStore, combineReducers } from '@reduxjs/toolkit';

export const getInitialState = () => ({
  visualize: {
    source: 'public',
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
    display: 'within',
    withinErr: false,
    refErr: false,
    pubErr: false,
    debugR: [],
    withinSubmitOverlay: false,
    refSubmitOverlay: false,
    pubSubmitOverlay: false,
  },
  mutationalPattern: {
    proportion: '0.8',
    pattern: 'NCG>NTG',
    txtPath: '',
    plotPath: '',
    plotURL: '',
    barPath: '',
    barURL: '',
    display: true,
    err: false,
    debugR: [],
    loading: false,
  },
  profileComparison: {
    withinProfileType: '',
    withinSampleName1: '',
    withinSampleName2: '',
    sampleOptions: '',
    refProfileType: '',
    refSampleName: '',
    refSignatureSet: '',
    refSignatureSetOptions: [],
    refSignatures: [],
    filterSignatures: [],
    refCompare: '',
    searchFilter: '',
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
    display: 'within',
    refErr: false,
    pubErr: false,
    debugR: [],
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
    pcaErr: false,
    debugR: [],
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
    display: 'within',
    pubPcaErr: false,
    pubSubmitOverlay: false,
  },
  kataegis: {
    sample: '',
    sampleOptions: '',
    highlight: false,
    min: '5',
    max: '100',
    chromosome: 'None',
    txtPath: '',
    plotPath: '',
    plotURL: '',
    kataegisData: [],
    dataField: '',
    order: '',
    page: 1,
    size: 10,
    display: true,
    err: false,
    debugR: [],
    loading: false,
  },

  exploring: {
    displayTab: 'etiology',
    signatureDisplay: 'referenceSignatures',
    exposureSignature: [],
    exposureCancer: [],
    refSigData: [],
    projectID: '',
    submitted: false,
  },
  expRefSig: {
    plotPath: '',
    plotURL: 'assets/images/mSigPortalReferenceSignatures.svg',
    debugR: '',
    err: '',
    loading: false,
  },
  expMutationalProfiles: {
    plots: [
      {
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
      },
    ],
    debugR: '',
    err: '',
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
    loading: false,
  },
  expExposure: {
    study: '',
    studyOptions: [],
    strategy: '',
    strategyOptions: [],
    cancer: '',
    cancerOptions: [],
    refSigData: [],
    refSignatureSet: '',
    refSignatureSetOptions: [],
    signatureNameOptions: [],
    userNameOptions: [],
    genome: 'GRCh37',
    genomeOptions: ['GRCh37', 'GRCh38', 'mm10'],
    exposureFile: '',
    matrixFile: '',
    signatureFile: '',
    usePublicSignature: true,
    source: 'public',
    display: 'tmb',
    loading: false,
    loadingMsg: null,
    projectID: '',
    openSidebar: true,
  },
  tmb: {
    plotPath: '',
    plotURL: '',
    debugR: '',
    err: '',
  },
  tmbSignatures: {
    plotPath: '',
    plotURL: '',
    debugR: '',
    err: '',
  },
  msBurden: {
    signatureName: '',
    plotPath: '',
    plotURL: '',
    debugR: '',
    err: '',
    loading: false,
  },
  msAssociation: {
    toggleCancer: true,
    both: true,
    signatureName1: '',
    signatureName2: '',
    plotPath: '',
    plotURL: '',
    txtPath: '',
    debugR: '',
    err: '',
    loading: false,
  },
  msDecomposition: {
    plotPath: '',
    plotURL: '',
    txtPath: '',
    debugR: '',
    err: '',
    loading: false,
  },
  msLandscape: {
    variableFile: '',
    plotPath: '',
    plotURL: '',
    txtPath: '',
    debugR: '',
    err: '',
    loading: false,
  },
  msPrevalence: {
    mutation: '100',
    plotPath: '',
    plotURL: '',
    debugR: '',
    err: '',
    loading: false,
  },
  error: {
    visible: false,
    message: `An error occured when requesting data. If this problem persists, please contact the administrator at <a href="mailto:mSigPortalWebAdmin@cancer.gov">mSigPortalWebAdmin@cancer.gov</a>.`,
  },
  success: {
    visible: false,
    message: 'Your job was successfuly submitted to the queue.',
  },
});

// visualization modules
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

const kataegisSlice = createSlice({
  name: 'kataegis',
  initialState: getInitialState().kataegis,
  reducers: {
    updateKataegis: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});

// util
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

const successSlice = createSlice({
  name: 'success',
  initialState: getInitialState().success,
  reducers: {
    updateSuccess: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});

// exploring modules
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
const expExposureSlice = createSlice({
  name: 'expExposure',
  initialState: getInitialState().tmb,
  reducers: {
    updateExpExposure: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});
const tmbSlice = createSlice({
  name: 'tmb',
  initialState: getInitialState().tmb,
  reducers: {
    updateTMB: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});
const tmbSignaturesSlice = createSlice({
  name: 'tmbSignatures',
  initialState: getInitialState().tmbSignatures,
  reducers: {
    updateTmbSignatures: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});
const msBurdenSlice = createSlice({
  name: 'msBurden',
  initialState: getInitialState().msBurden,
  reducers: {
    updateMsBurden: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});
const msAssocationSlice = createSlice({
  name: 'msAssocation',
  initialState: getInitialState().msAssociation,
  reducers: {
    updateMsAssociation: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});
const msDecompositionSlice = createSlice({
  name: 'msDecomposition',
  initialState: getInitialState().msDecomposition,
  reducers: {
    updateMsDecomposition: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});
const msLandscapeSlice = createSlice({
  name: 'msLandscape',
  initialState: getInitialState().msLandscape,
  reducers: {
    updateMsLandscape: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});
const msPrevalenceSlice = createSlice({
  name: 'msPrevalence',
  initialState: getInitialState().msPrevalence,
  reducers: {
    updateMsPrevalence: (state, action) => {
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
  kataegis: kataegisSlice.reducer,

  error: errorSlice.reducer,
  success: successSlice.reducer,

  exploring: exploringSlice.reducer,
  expRefSig: expRefSigSlice.reducer,
  expMutationalProfiles: expMutationalProfilesSlice.reducer,
  expCosineSimilarity: expCosineSimilaritySlice.reducer,
  expMutationalSigComparison: expMutationalSigComparisonSlice.reducer,
  expExposure: expExposureSlice.reducer,
  tmb: tmbSlice.reducer,
  tmbSignatures: tmbSignaturesSlice.reducer,
  msBurden: msBurdenSlice.reducer,
  msAssociation: msAssocationSlice.reducer,
  msDecomposition: msDecompositionSlice.reducer,
  msLandscape: msLandscapeSlice.reducer,
  msPrevalence: msPrevalenceSlice.reducer,
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
export const { updateKataegis } = kataegisSlice.actions;

export const { updateError } = errorSlice.actions;
export const { updateSuccess } = successSlice.actions;

export const { updateExploring } = exploringSlice.actions;
export const { updateExpRefSig } = expRefSigSlice.actions;
export const {
  updateExpMutationalProfiles,
} = expMutationalProfilesSlice.actions;
export const { updateExpCosineSimilarity } = expCosineSimilaritySlice.actions;
export const {
  updateExpMutationalSigComparison,
} = expMutationalSigComparisonSlice.actions;
export const { updateExpExposure } = expExposureSlice.actions;
export const { updateTMB } = tmbSlice.actions;
export const { updateTmbSignatures } = tmbSignaturesSlice.actions;
export const { updateMsBurden } = msBurdenSlice.actions;
export const { updateMsAssociation } = msAssocationSlice.actions;
export const { updateMsDecomposition } = msDecompositionSlice.actions;
export const { updateMsLandscape } = msLandscapeSlice.actions;
export const { updateMsPrevalence } = msPrevalenceSlice.actions;

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

export function dispatchKataegis(obj) {
  store.dispatch(updateKataegis(obj));
}

export function dispatchError(msg) {
  store.dispatch(
    updateError({
      visible: true,
      message: msg,
    })
  );
}

export function dispatchSuccess(msg) {
  store.dispatch(
    updateSuccess({
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
export function dispatchExpExposure(obj) {
  store.dispatch(updateExpExposure(obj));
}
export function dispatchTMB(obj) {
  store.dispatch(updateTMB(obj));
}
export function dispatchTmbSignatures(obj) {
  store.dispatch(updateTmbSignatures(obj));
}
export function dispatchMsBurden(obj) {
  store.dispatch(updateMsBurden(obj));
}
export function dispatchMsAssociation(obj) {
  store.dispatch(updateMsAssociation(obj));
}
export function dispatchMsDecomposition(obj) {
  store.dispatch(updateMsDecomposition(obj));
}
export function dispatchMsLandscape(obj) {
  store.dispatch(updateMsLandscape(obj));
}
export function dispatchMsPrevalence(obj) {
  store.dispatch(updateMsPrevalence(obj));
}
