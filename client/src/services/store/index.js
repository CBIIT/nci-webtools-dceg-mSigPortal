import { reducer as catalogReducer } from './catalog';
import { reducer as exposureReducer } from './exposure';
import { reducer as visualizationReducer } from './visualization';
import { reducer as associationReducer } from './association';
import { reducer as modalReducer } from './modal';
import { reducer as publicationsReducer } from './publications';

import { configureStore } from '@reduxjs/toolkit';
import { setupListeners } from '@reduxjs/toolkit/query';
import {
  optionsApiSlice,
  visualizationApiSlice,
  explorationApiSlice,
} from './rootApi';

// provide rootReducer as an object of slice reducers
export const store = configureStore({
  reducer: {
    visualization: visualizationReducer,
    catalog: catalogReducer,
    association: associationReducer,
    exposure: exposureReducer,
    modal: modalReducer,
    publications: publicationsReducer,

    [optionsApiSlice.reducerPath]: optionsApiSlice.reducer,
    [visualizationApiSlice.reducerPath]: visualizationApiSlice.reducer,
    [explorationApiSlice.reducerPath]: explorationApiSlice.reducer,
  },
  middleware: (getDefaultMiddleware) =>
    getDefaultMiddleware({ serializableCheck: false }).concat(
      optionsApiSlice.middleware,
      visualizationApiSlice.middleware,
      explorationApiSlice.middleware
    ),
});

setupListeners(store.dispatch);
