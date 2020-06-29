import { createSlice,configureStore} from '@reduxjs/toolkit';
// import { getInitialState } from './store-old'

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

const visualSlice = createSlice({
    name: 'visualize',
    initialState: getInitialState().visualize,
    reducers: {
        setInputFormat: (state,action) => {
            return{
                ...state,
                inputFormat: action.payload
            }
        } 
    }
  })

const {actions,reducer} = visualSlice
export const store = configureStore({
    reducer: reducer,
    preloadedState: getInitialState()
})
export const {setInputFormat} = actions
// export default store
