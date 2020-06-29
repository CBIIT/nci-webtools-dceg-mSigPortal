import { createSlice,configureStore} from '@reduxjs/toolkit';

const counterSlice = createSlice({

    name:'counter',
    initialState: 0,
    reducers:{
        increment: state => state + 1,
        decrement: state => state - 1
      }
})

const {actions,reducer} = counterSlice
const store = configureStore({reducer: reducer})

export const {increment,decrement} = actions
export default store
