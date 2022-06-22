import { createSlice } from '@reduxjs/toolkit';
import { mergeObject } from './utils';

export const getInitialState = () => ({
  publicForm: {
    defaultOptions: {
      study: { label: 'PCAWG', value: 'PCAWG' },
      cancer: { label: 'Lung-AdenoCA', value: 'Lung-AdenoCA' },
      strategy: { label: 'WGS', value: 'WGS' },
    },
    study: null,
    cancer: null,
    strategy: null,
    loading: false,
    data: [],
  },
  userForm: {
    inputFormat: 'vcf',
    selectedGenome: 'GRCh37',
    experimentalStrategy: 'WGS',
    mutationSplit: 'False',
    collapseSample: 'False',
    mutationFilter: '',
    queueMode: false,
    email: '',
  },
  main: {
    submitted: false,
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
    exampleData: 'assets/exampleInput/demo_input_multi.vcf.gz',
    bedData: 'assets/exampleInput/demo_input_bed.bed',
    loadingPublic: false,
    loading: {
      active: false,
      content: null,
      showIndicator: false,
    },
    queueExpired: false,
    error: '',
    projectID: '',
    displayTab: 'instructions',
    downloads: {},
    svgList: {},
    matrixList: {},
    statistics: '',
    profileOptions: [],
  },
  profilerSummary: {
    plotPath: '',
    err: '',
    debugR: '',
    loading: false,
  },
  mutationalProfiles: {
    sample: null,
    profile: null,
    matrix: null,
    filter: null,
    data: [],
    plot: null,
    loading: false,

    filtered: [],
    selectName: '',
    selectProfile: '',
    selectMatrix: '',
    selectFilter: '',
    nameOptions: [],
    profileOptions: [],
    matrixOptions: [],
    filterOptions: [],
    plotPath: '',
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
    proportion: '',
    pattern: '',
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
    pubPlotPath: '',
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
    kataegisData: [],
    pagination: {
      pageIndex: 0,
      pageSize: 10,
    },
    hidden: [],
    display: true,
    err: false,
    debugR: [],
    loading: false,
  },
});

export const { actions, reducer } = createSlice({
  name: 'visualization',
  initialState: getInitialState(),
  reducers: {
    mergeVisualization: mergeObject,
    resetVisualization: getInitialState,
  },
});
