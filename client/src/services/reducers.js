import {
  UPDATE_VISUALIZE,
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
    default:
      return state;
  }
};
