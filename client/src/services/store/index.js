import { reducer as explorationReducer } from './exploration';
import { reducer as visualizationReducer } from './visualization';
import { reducer as associationReducer } from './association';
import { reducer as modalReducer } from './modal';
import { reducer as publicationsReducer } from './publications';

import { configureStore, getDefaultMiddleware } from '@reduxjs/toolkit';

// provide rootReducer as an object of slice reducers
export const store = configureStore({
  reducer: {
    visualization: visualizationReducer,
    exploration: explorationReducer,
    association: associationReducer,
    modal: modalReducer,
    publications: publicationsReducer,
  },
  middleware: getDefaultMiddleware({
    serializableCheck: false,
  }),
});
