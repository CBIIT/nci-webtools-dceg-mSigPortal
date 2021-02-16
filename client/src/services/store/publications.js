import { createSlice } from '@reduxjs/toolkit';
import { mergeObject } from './utils';

export const getInitialState = () => ({});

export const { actions, reducer } = createSlice({
  name: 'publications',
  initialState: getInitialState(),
  reducers: {
    mergeState: mergeObject,
    resetModal: getInitialState,
  },
});
