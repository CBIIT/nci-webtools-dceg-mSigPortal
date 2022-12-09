import { createSlice } from '@reduxjs/toolkit';
import { mergeObject } from './utils';

export const getInitialState = () => ({
  main: {
    displayTab: 'etiology',

    signatureDisplay: 'referenceSignatures',
    exposureSignature: [],
    exposureCancer: [],
    refSigData: [],
    projectID: '',
    submitted: false,
  },
  etiology: {
    category: 'Cosmic',
    etiology: '',
    signature: '',
    referrenceSignature: '',
    study: '',
    cohort: '',
    showAllSignatures: false,
  },
  referenceSignature: {
    displayTab: 'overview',
    refSigData: [],
    exposureSignature: [],
    projectID: '',
  },
  sigRefSig: {
    plot: '',
    err: '',
  },
  rSProfiles: {
    plots: [
      {
        source: '',
        profile: '',
        matrix: '',
        signatureSetName: '',
        strategy: '',
        signatureName: '',
      },
    ],
  },
  cosineSimilarity: {
    profile: '',
    matrix: '',
    signatureSet1: '',
    signatureSet2: '',
  },
  sigMutationalSigComparison: {
    profileName: '',
    profileNameOptions: [],
    rsSet1: '',
    rsSetOptions1: [],
    rsSet2: '',
    rsSetOptions2: [],
    signatureName1: '',
    signatureNameOptions1: [],
    signatureName2: '',
    signatureNameOptions2: [],
    plotPath: '',
    debugR: '',
    err: '',
    loading: false,
  },
});

export const { actions, reducer } = createSlice({
  name: 'catalog',
  initialState: getInitialState(),
  reducers: {
    mergeCatalog: mergeObject,
    resetCatalog: getInitialState,
  },
});
