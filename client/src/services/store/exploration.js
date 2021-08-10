import { createSlice } from '@reduxjs/toolkit';
import { mergeObject } from './utils';

export const getInitialState = () => ({
  exploration: {
    displayTab: 'etiology',
    signatureDisplay: 'referenceSignatures',
    exposureSignature: [],
    exposureCancer: [],
    refSigData: [],
    projectID: '',
    submitted: false,
  },
  sigRefSig: {
    plotPath: 'assets/images/mSigPortalReferenceSignatures.svg',
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
    debugR: '',
    err: '',
    loading: false,
  },
});

export const { actions, reducer } = createSlice({
  name: 'exploration',
  initialState: getInitialState(),
  reducers: {
    mergeExploration: mergeObject,
    resetExploration: getInitialState,
  },
});
