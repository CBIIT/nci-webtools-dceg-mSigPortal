import { createSlice, configureStore, combineReducers } from '@reduxjs/toolkit';

export const getInitialState = () => ({
    visualize: {
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
    },
    error: {
        visible: false,
        message: `An error occured when requesting data. If this problem persists, please contact the administrator at <a href="mailto:mSigPortalWebAdmin@cancer.gov">mSigPortalWebAdmin@cancer.gov</a>.`
    }
});

const visualizeSlice = createSlice({
    name: 'visualize',
    initialState: getInitialState().visualize,
    reducers: {
        replaceVisualize: (state, action) => {
            return {
                ...state,
                ...action.payload.data
            }
        },
        updateVisualize: (state, action) => {
            return {
                ...state,
                [action.payload.param]: action.payload.data
            }
        },
    }
});

const visualizeResultsSlice = createSlice({
    name: 'visualizeResults',
    initialState: getInitialState().visualizeResults,
    reducers: {
        replaceVisualizeResults: (state, action) => {
            return {
                ...state,
                ...action.payload.data
            }
        },
        updateVisualizeResults: (state, action) => {
            return {
                ...state,
                [action.payload.param]: action.payload.data
            }
        }
    }
});

const errorSlice = createSlice({
    name: 'error',
    initialState: getInitialState().error,
    reducers: {
        updateError: (state, action) => {
            return {
                ...state,
                [action.payload.param]: action.payload.data
            }
        }
    }
});

const rootReducer = combineReducers({
    visualize: visualizeSlice.reducer,
    visualizeResults: visualizeResultsSlice.reducer,
    error: errorSlice.reducer
})

export const store = configureStore({
    reducer: rootReducer,
    preloadedState: getInitialState()
});

export const {
    replaceVisualize,
    updateVisualize, 
} = visualizeSlice.actions;

export const {
    replaceVisualizeResults,
    updateVisualizeResults
} = visualizeResultsSlice.actions;

export const {
    updateError
} = errorSlice.actions;
