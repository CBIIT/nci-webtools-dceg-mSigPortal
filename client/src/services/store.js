import { rootReducer } from './reducers';
// import { initialize } from './actions';
import {
  // createAction,
  configureStore,
} from '@reduxjs/toolkit';

export const getInitialState = () => ({
  visualizeForm: {
    count: 0,
    inputFormat: 'vcf',
    selectedGenome: 'GRCh37',
    experimentalStrategy: 'WGS',
    mutationSplit: 'False',
    isMultiple: false,
    collapseSample: 'False',
    mutationFilter: '',
    queueMode: false,
    email: '',
  },
  visualizeResults: {
    projectID: '',
    mapping: [],
    plots: [],
    displayedPlot: '',
    plotURL: '',
  },
});

export const store = configureStore({
  reducer: rootReducer,
  preloadedState: getInitialState(),
});

// store.dispatch(initialize());
