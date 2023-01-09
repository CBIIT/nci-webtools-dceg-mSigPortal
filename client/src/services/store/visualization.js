import { createSlice } from '@reduxjs/toolkit';
import { mergeObject } from './utils';

export const getInitialState = () => ({
  publicForm: {
    study: '',
    cancer: '',
    strategy: '',
  },
  userForm: {
    inputFilename: '',
    bedFilename: '',
    inputFormat: {
      label: 'VCF',
      value: 'vcf',
      example: 'demo_input_multi.vcf.gz',
    },
    genome: '',
    strategy: '',
    mutationSplit: '',
    collapse: false,
    filter: '',
    useQueue: false,
    email: '',
  },
  main: {
    matrixData: [],
    submitted: false,
    source: 'public',
    openSidebar: true,

    loading: {
      active: false,
      content: '',
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
  },
  profilerSummary: { filter: '' },
  mutationalProfiles: {
    sample: '',
    profile: '',
    matrix: '',
    filter: '',
    data: [],
    plot: '',
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
    withinForm: { profile: '', matrix: '' },
    referenceForm: { profile: '', signatureSet: '' },
    publicForm: {},

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
    withinForm: { profile: '', sample1: '', sample2: '' },
    referenceForm: {
      profile: { value: 'SBS', label: 'SBS' },
      sample: { value: 'SP50263', label: 'SP50263' },
      signatureSet: {
        value: 'COSMIC_v3.3_Signatures_GRCh38_SBS96',
        label: 'COSMIC_v3.3_Signatures_GRCh38_SBS96',
      },
      compare: { value: 'SBS1', label: 'SBS1' },
    },
    publicForm: {},
    display: 'within',
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
  clustered: {
    sample: '',
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
