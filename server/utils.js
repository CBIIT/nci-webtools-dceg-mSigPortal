// get row value from index of column key
function value2d(row, key, columns) {
  return row[columns.indexOf(key)];
}

// filter 2d array by a given value
function filter2d(value, data) {
  return Array.isArray(value)
    ? data.filter((row) => value.every((v) => row.includes(v)))
    : data.filter((row) => row.includes(value));
}

// get unique array of row entries for a given key (column) from 2d array
function unique2d(key, columns, data) {
  return [...new Set(data.map((row) => row[columns.indexOf(key)]))].sort(
    (a, b) => a - b
  );
}

function defaultProfile(profileOptions) {
  if (profileOptions.includes('SBS')) return 'SBS';
  if (profileOptions.includes('DBS')) return 'DBS';
  if (profileOptions.includes('ID')) return 'ID';
}

function defaultMatrix(profile, matrixOptions) {
  if (profile == 'SBS')
    return matrixOptions.includes('96') ? '96' : matrixOptions[0];

  if (profile == 'DBS')
    return matrixOptions.includes('78') ? '78' : matrixOptions[0];

  if (profile == 'ID')
    return matrixOptions.includes('83') ? '83' : matrixOptions[0];
}

function defaultFilter(filterOptions) {
  return filterOptions.includes('NA') ? 'NA' : filterOptions[0];
}

module.exports = {
  value2d,
  filter2d,
  unique2d,
  defaultProfile,
  defaultMatrix,
  defaultFilter,
};
