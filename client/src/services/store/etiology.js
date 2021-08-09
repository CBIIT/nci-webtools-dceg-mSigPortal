import { createSlice } from '@reduxjs/toolkit';
import { mergeObject } from './utils';

export const getInitialState = () => ({
  etiologyState: {
    category: 'Cosmic Mutational Signatures',
    etiology: 'APOBEC activity',
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
});

export const { actions, reducer } = createSlice({
  name: 'etiology',
  initialState: getInitialState(),
  reducers: {
    mergeEtiology: mergeObject,
    resetEtiology: getInitialState,
  },
});
