import { rootReducer } from './reducers';
// import { initialize } from './actions';
import { 
  // createAction,
  configureStore
} from '@reduxjs/toolkit';

export const getInitialState = () => ({
  visualize: {
    count: 0,
    inputFormat: 'vcf',
    inputFile: null,
    selectedGenome: 'GRCh37',
    experimentalStrategy: 'WGS',
    mutationSplit: 'False',
    isMultiple: false,
    collapseSample: 'False',
    mutationFilter: '',
    queueMode: false,
    email: ''
  },
  visualizeResults: {
    uid: null,
    mapping: null,
    displayedPlot: null
  }
});

export const store = configureStore({
  reducer: rootReducer,
  preloadedState: getInitialState()
});

// store.dispatch(initialize());
