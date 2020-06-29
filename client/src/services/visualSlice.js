import { createSlice,configureStore} from '@reduxjs/toolkit';
import { getInitialState } from './store'

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
const store = configureStore({reducer: reducer})
export const {setInputFormat} = actions
export default store
