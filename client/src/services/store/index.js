import { reducer as catalogReducer } from './catalog';
import { reducer as visualizationReducer } from './visualization';
import { reducer as explorationReducer } from './exploration';
import { reducer as associationReducer } from './association';
import { reducer as extractionReducer } from './extraction';
import { reducer as refittingReducer } from './refitting';
import { reducer as modalReducer } from './modal';

import { configureStore } from '@reduxjs/toolkit';
import { setupListeners } from '@reduxjs/toolkit/query';
import {
  optionsApiSlice,
  visualizationApiSlice,
  explorationApiSlice,
  catalogApiSlice,
  associationApiSlice,
  extractionApiSlice,
  refittingApiSlice,
} from './rootApi';

// provide rootReducer as an object of slice reducers
export const store = configureStore({
  reducer: {
    catalog: catalogReducer,
    visualization: visualizationReducer,
    exploration: explorationReducer,
    association: associationReducer,
    extraction: extractionReducer,
    refitting: refittingReducer,
    modal: modalReducer,

    [optionsApiSlice.reducerPath]: optionsApiSlice.reducer,
    [visualizationApiSlice.reducerPath]: visualizationApiSlice.reducer,
    [explorationApiSlice.reducerPath]: explorationApiSlice.reducer,
    [catalogApiSlice.reducerPath]: catalogApiSlice.reducer,
    [associationApiSlice.reducerPath]: associationApiSlice.reducer,
    [extractionApiSlice.reducerPath]: extractionApiSlice.reducer,
    [refittingApiSlice.reducerPath]: refittingApiSlice.reducer,
  },
  middleware: (getDefaultMiddleware) =>
    getDefaultMiddleware({ serializableCheck: false }).concat(
      optionsApiSlice.middleware,
      visualizationApiSlice.middleware,
      explorationApiSlice.middleware,
      catalogApiSlice.middleware,
      associationApiSlice.middleware,
      extractionApiSlice.middleware,
      refittingApiSlice.middleware
    ),
});

setupListeners(store.dispatch);
