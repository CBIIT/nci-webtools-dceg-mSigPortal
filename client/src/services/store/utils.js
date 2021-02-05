import _ from 'lodash';

export function mergeObject(state, action) {
  return _.merge(state, action.payload);
}

export function mergeArray(state, action) {
  return [...state, ...action.payload];
}

export function pushArray(state, action) {
  return [...state, action.payload];
}
