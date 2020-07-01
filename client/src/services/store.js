import { createSlice, configureStore, combineReducers } from '@reduxjs/toolkit';

export const getInitialState = () => ({
  visualize: {
    inputFormat: 'vcf',
    selectedGenome: 'GRCh37',
    experimentalStrategy: 'WGS',
    mutationSplit: 'False',
    isMultiple: false,
    collapseSample: 'False',
    mutationFilter: '',
    queueMode: false,
    email: '',
    openSidebar: true,
    disableParameters: false,
    loading: {
      active: false,
      content: null,
      showIndicator: false
    }
  },
  visualizeResults: {
    error: '',
    projectID: '',
    mapping: [],
    plots: [],
    displayedPlot: '',
    plotURL: '',
  },
  error: {
    visible: false,
    message: `An error occured when requesting data. If this problem persists, please contact the administrator at <a href="mailto:mSigPortalWebAdmin@cancer.gov">mSigPortalWebAdmin@cancer.gov</a>.`,
  },
});

const visualizeSlice = createSlice({
  name: 'visualize',
  initialState: getInitialState().visualize,
  reducers: {
    updateVisualize: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});

const visualizeResultsSlice = createSlice({
  name: 'visualizeResults',
  initialState: getInitialState().visualizeResults,
  reducers: {
    updateVisualizeResults: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});

const errorSlice = createSlice({
  name: 'error',
  initialState: getInitialState().error,
  reducers: {
    updateError: (state, action) => {
      return {
        ...state,
        ...action.payload,
      };
    },
  },
});

const rootReducer = combineReducers({
  visualize: visualizeSlice.reducer,
  visualizeResults: visualizeResultsSlice.reducer,
  error: errorSlice.reducer,
});

export const store = configureStore({
  reducer: rootReducer,
  preloadedState: getInitialState(),
});

export const { updateVisualize } = visualizeSlice.actions;

export const { updateVisualizeResults } = visualizeResultsSlice.actions;

export const { updateError } = errorSlice.actions;
