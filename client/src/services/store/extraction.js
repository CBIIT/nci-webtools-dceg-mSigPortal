import { createSlice } from '@reduxjs/toolkit';
import { mergeObject } from './utils';

export const getInitialState = () => ({
  tabIndex: 0,
  id: '',
  openSidebar: true,
  submitted: false,
});

export const { actions, reducer } = createSlice({
  name: 'extraction',
  initialState: getInitialState(),
  reducers: {
    mergeExtraction: mergeObject,
    resetExtraction: getInitialState,
  },
});
