import { createSlice,configureStore} from '@reduxjs/toolkit';
// import { getInitialState } from './store-old'

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
    }
  });

const visualSlice = createSlice({
    name: 'visualize',
    initialState: getInitialState(),
    reducers: {
        updateVisualize: (state,action) => {
            return{
                ...state,
                visualize: {
                    ...state.visualize,
                    [action.payload.param]: action.payload.data
                }
            }
        },
        resetVisualize: (state) => {
            return{
                ...state,
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
                }
            }
        }
    }
  })

const {actions,reducer} = visualSlice
export const store = configureStore({
    reducer: reducer,
    preloadedState: getInitialState()
})
export const {updateVisualize,resetVisualize} = actions
// export default store
