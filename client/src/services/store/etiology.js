import { createSlice } from '@reduxjs/toolkit';
import { mergeObject } from './utils';

export const getInitialState = () => ({
  category: 'Cosmic Mutational Signatures (v3.2)',
  etiology: 'APOBEC activity',
  signatureName: '',
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
  strandbiasURL: '',
  tissueURL: '',
  refSigURL: '',
});

export const { actions, reducer } = createSlice({
  name: 'etiology',
  initialState: getInitialState(),
  reducers: {
    mergeEtiology: mergeObject,
    resetEtiology: getInitialState,
  },
});
