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
  },
  visualizeResults: {
    projectID: '',
    mapping: [],
    plots: [],
    displayedPlot: '',
    plotURL: '',
  },
});

const visualSlice = createSlice({
  name: 'visualize',
  initialState: getInitialState(),
  reducers: {
    replaceVisualize: (state, action) => {
      return {
        ...state,
        [action.payload.param]: action.payload.data,
      };
    },
    updateVisualize: (state, action) => {
      return {
        ...state,
        visualize: {
          ...state.visualize,
          [action.payload.param]: action.payload.data,
        },
      };
    },
    updateVisualizeResults: (state, action) => {
      return {
        ...state,
        visualizeResults: {
          ...state.visualizeResults,
          [action.payload.param]: action.payload.data,
        },
      };
    },
    // resetVisualize: (state) => {
    //     return{
    //         ...state,
    //         visualize: {
    //             inputFormat: 'vcf',
    //             inputFile: null,
    //             selectedGenome: 'GRCh37',
    //             experimentalStrategy: 'WGS',
    //             mutationSplit: 'False',
    //             isMultiple: false,
    //             collapseSample: 'False',
    //             mutationFilter: '',
    //             queueMode: false,
    //             email: ''
    //         }
    //     }
    // }
  },
});

const { actions, reducer } = visualSlice;
export const store = configureStore({
  reducer: reducer,
  preloadedState: getInitialState(),
});

export const {
  replaceVisualize,
  updateVisualize,
  updateVisualizeResults,
  // resetVisualize
} = actions;
