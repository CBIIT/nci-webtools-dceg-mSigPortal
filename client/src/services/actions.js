// import { 
//   createAction,
// } from '@reduxjs/toolkit';
// import { getInitialState } from './store';

export const UPDATE_VISUALIZE = 'UPDATE_VISUALIZE';

export function updateVisualize(data) {
  return { type: UPDATE_VISUALIZE, data };
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


