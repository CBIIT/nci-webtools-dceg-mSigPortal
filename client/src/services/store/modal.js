import { createSlice } from '@reduxjs/toolkit';
import { mergeObject } from './utils';

export const getInitialState = () => ({
  error: {
    visible: false,
    message: [],
  },
  success: {
    visible: false,
    message: 'Your job was successfuly submitted to the queue.',
  },
});

export const { actions, reducer } = createSlice({
  name: 'modal',
  initialState: getInitialState(),
  reducers: {
    mergeModal: mergeObject,
    resetModal: getInitialState,
  },
});
