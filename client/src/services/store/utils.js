import _ from 'lodash';

export function mergeObject(state, action) {
  return _.mergeWith(state, action.payload, (obj, src) => {
    if (_.isArray(obj)) {
      return src;
    }
  });
}

export function mergeArray(state, action) {
  return [...state, ...action.payload];
}

export function pushArray(state, action) {
  return [...state, action.payload];
}
