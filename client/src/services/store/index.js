import { reducer as exploringReducer } from './exploring';
import { reducer as visualizationReducer } from './visualization';
import { reducer as modalReducer } from './modal';
import { reducer as publicationsReducer } from './publications';

import { configureStore, getDefaultMiddleware } from '@reduxjs/toolkit';

// provide rootReducer as an object of slice reducers
export const store = configureStore({
  reducer: {
    visualization: visualizationReducer,
    exploring: exploringReducer,
    modal: modalReducer,
    publications: publicationsReducer,
  },
  middleware: getDefaultMiddleware({
    serializableCheck: false,
  }),
});
