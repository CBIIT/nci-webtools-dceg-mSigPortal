import { reducer as exploringReducer } from './exploring';
import { reducer as visualizationReducer } from './visualization';
import { reducer as modalReducer } from './modal';

import { configureStore } from '@reduxjs/toolkit';

// provide rootReducer as an object of slice reducers
export const store = configureStore({
  reducer: {
    visualization: visualizationReducer,
    exploring: exploringReducer,
    modal: modalReducer,
  },
});
