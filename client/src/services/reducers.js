import {
  UPDATE_VISUALIZE,
  UPDATE_VISUALIZE_RESULTS
} from './actions';

export const rootReducer = (state, action) => {
  switch (action.type) {
    case UPDATE_VISUALIZE:
      return {
        ...state,
        visualize: {
          ...state.visualize,
          ...action.data
        }
      };
    case UPDATE_VISUALIZE_RESULTS:
      return {
        ...state,
        visualizeResults: {
          ...state.visualizeResults,
          ...action.data
        }
      };
    default:
      return state;
  }
};
