import { createSlice } from '@reduxjs/toolkit';
import { mergeObject } from './utils';

export const getInitialState = () => ({
  error: {
    visible: false,
    message: `An error occured when requesting data. If this problem persists, please contact the administrator at <a href="mailto:mSigPortalWebAdmin@cancer.gov">mSigPortalWebAdmin@cancer.gov</a>.`,
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
