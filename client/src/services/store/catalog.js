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
        profile: { value: 'SBS', label: 'SBS' },
        matrix: { value: '96', label: '96' },
        signatureSetName: {
          value: 'Other_Published_Signatures_GRCh37_SBS96',
          label: 'Other_Published_Signatures_GRCh37_SBS96',
        },
        strategy: { value: 'TS', label: 'TS' },
        signatureName: { value: 'SBS_CPD_CDSEQ', label: 'SBS_CPD_CDSEQ' },
      },
    ],
  },
  cosineSimilarity: {
    profile: { value: 'SBS', label: 'SBS' },
    matrix: { value: '96', label: '96' },
    signatureSet1: {
      value: 'COSMIC_v3.3_Signatures_GRCh38_SBS96',
      label: 'COSMIC_v3.3_Signatures_GRCh38_SBS96',
    },
    signatureName1: {
      value: 'SBS1',
      label: 'SBS1',
    },
    signatureSet2: {
      value: 'COSMIC_v1_Signatures_GRCh37_SBS96',
      label: 'COSMIC_v1_Signatures_GRCh37_SBS96',
    },
    signatureName2: {
      value: 'Signature_1A',
      label: 'Signature_1A',
    },
  },
  sigMutationalSigComparison: {
    profileName: { value: 'SBS', label: 'SBS' },
    profileNameOptions: [],
    rsSet1: {
      value: 'COSMIC_v3.3_Signatures_GRCh38_SBS96',
      label: 'COSMIC_v3.3_Signatures_GRCh38_SBS96',
    },
    rsSetOptions1: [],
    rsSet2: {
      value: 'COSMIC_v1_Signatures_GRCh37_SBS96',
      label: 'COSMIC_v1_Signatures_GRCh37_SBS96',
    },
    rsSetOptions2: [],
    signatureName1: {
      value: 'SBS1',
      label: 'SBS1',
    },
    signatureNameOptions1: [],
    signatureName2: {
      value: 'Signature_1A',
      label: 'Signature_1A',
    },
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
