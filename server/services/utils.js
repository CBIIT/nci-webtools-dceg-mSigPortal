import { pickBy } from 'lodash-es';

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
function defaultProfile2(profileOptions) {
  const options = profileOptions.map(({ value }) => value);
  if (options.includes('SBS')) return { label: 'SBS', value: 'SBS' };
  if (options.includes('DBS')) return { label: 'DBS', value: 'DBS' };
  if (options.includes('ID')) return { label: 'ID', value: 'ID' };
}
function defaultMatrix2(profile, matrixOptions) {
  const options = options.map(({ value }) => value);
  if (profile.value == 'SBS')
    return options.includes('96')
      ? { label: '96', value: '96' }
      : matrixOptions[0];
  if (profile.value == 'DBS')
    return options.includes('78')
      ? { label: '78', value: '78' }
      : matrixOptions[0];
  if (profile.value == 'ID')
    return options.includes('83')
      ? { label: '83', value: '83' }
      : matrixOptions[0];
}
function defaultFilter2(filterOptions) {
  const options = options.map(({ value }) => value);
  return options.includes('NA')
    ? { label: 'NA', value: 'NA' }
    : filterOptions[0];
}
function pickNonNullValues(object) {
  return pickBy(object, (v) => v !== null);
}

export {
  defaultProfile,
  defaultMatrix,
  defaultFilter,
  defaultProfile2,
  defaultMatrix2,
  defaultFilter2,
  pickNonNullValues,
};
