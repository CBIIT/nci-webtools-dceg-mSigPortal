import { createSlice } from '@reduxjs/toolkit';
import { mergeObject } from './utils';

export const getInitialState = () => ({
  exploring: {
    displayTab: 'aetiology',
    signatureDisplay: 'referenceSignatures',
    exposureSignature: [],
    exposureCancer: [],
    refSigData: [],
    projectID: '',
    submitted: false,
  },
  aetiology: {
    category: 'Cosmic Mutational Signatures',
    aetiology: 'APOBEC activity',
    signature: '',
    tissue: '',
    refSig: '',
    study: 'PCAWG',
    all: false,
    data: [],
    thumbnails: [],
    tissueThumbnails: [],
    refSigThumbnails: [],
    selectedSource: '',
    profileURL: '',
    exposureURL: '',
    tissueURL: '',
    refSigURL: '',
  },
  sigRefSig: {
    plotPath: '',
    plotURL: 'assets/images/mSigPortalReferenceSignatures.svg',
    debugR: '',
    err: '',
    loading: false,
  },
  sigMutationalProfiles: {
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
  sigCosineSimilarity: {
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
  sigMutationalSigComparison: {
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
  exposure: {
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
    publicSampleOptions: [],
    userSampleOptions: [],
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
    loading: false,
  },
  tmbSignatures: {
    plotPath: '',
    plotURL: '',
    debugR: '',
    err: '',
    loading: false,
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
  msIndividual: {
    sample: '',
    plotPath: '',
    plotURL: '',
    debugR: '',
    err: '',
    loading: false,
  },
});

export const { actions, reducer } = createSlice({
  name: 'exploring',
  initialState: getInitialState(),
  reducers: {
    mergeExploring: mergeObject,
    resetExploring: getInitialState,
  },
});
