// import { 
//   createAction,
// } from '@reduxjs/toolkit';
// import { getInitialState } from './store';

export const UPDATE_VISUALIZE = 'UPDATE_VISUALIZE';
export const UPDATE_VISUALIZE_RESULTS = 'UPDATE_VISUALIZE_RESULTS';

export function updateVisualize(data) {
  return { type: UPDATE_VISUALIZE, data };
}

export function updateVisualizeResults(data) {
  return { type: UPDATE_VISUALIZE_RESULTS, data };
}

// const updateVisualize = createAction('UPDATE_VISUALIZE', function prepare(text) {
//   return {
//     payload: { type: UPDATE_VISUALIZE, data }
//   }
// })

// export function initialize() {
//   return async function initializeAction(dispatch) {
//     try {
//       console.log("Do something...");
//     } catch (e) {
// 
//     }
//   }
// }


