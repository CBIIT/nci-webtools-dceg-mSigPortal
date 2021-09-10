import { reducer as catalogReducer } from './catalog';
import { reducer as exposureReducer } from './exposure';
import { reducer as visualizationReducer } from './visualization';
import { reducer as associationReducer } from './association';
import { reducer as modalReducer } from './modal';
import { reducer as publicationsReducer } from './publications';

import { configureStore, getDefaultMiddleware } from '@reduxjs/toolkit';

// provide rootReducer as an object of slice reducers
export const store = configureStore({
  reducer: {
    visualization: visualizationReducer,
    catalog: catalogReducer,
    association: associationReducer,
    exposure: exposureReducer,
    modal: modalReducer,
    publications: publicationsReducer,
  },
  middleware: getDefaultMiddleware({
    serializableCheck: false,
  }),
});
