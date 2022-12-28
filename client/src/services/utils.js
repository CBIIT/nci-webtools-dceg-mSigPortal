export function defaultProfile(profileOptions) {
  if (profileOptions.includes('SBS')) return 'SBS';
  if (profileOptions.includes('DBS')) return 'DBS';
  if (profileOptions.includes('ID')) return 'ID';
}

export function defaultMatrix(profile, matrixOptions) {
  if (profile == 'SBS')
    return matrixOptions.includes('96') ? '96' : matrixOptions[0];

  if (profile == 'DBS')
    return matrixOptions.includes('78') ? '78' : matrixOptions[0];

  if (profile == 'ID')
    return matrixOptions.includes('83') ? '83' : matrixOptions[0];
}

export function defaultFilter(filterOptions) {
  return filterOptions.includes('NA') ? 'NA' : filterOptions[0];
}

export function defaultProfile2(profileOptions) {
  const options = profileOptions.map(({ value }) => value);
  if (options.includes('SBS')) return { label: 'SBS', value: 'SBS' };
  if (options.includes('DBS')) return { label: 'DBS', value: 'DBS' };
  if (options.includes('ID')) return { label: 'ID', value: 'ID' };
  if (options.includes('RS')) return { label: 'RS', value: 'RS' };
  if (options.includes('CN')) return { label: 'CN', value: 'CN' };
}

export function defaultMatrix2(profile, matrixOptions) {
  const options = matrixOptions.map(({ value }) => value);

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

  if (profile.value == 'RS')
    return options.includes('32')
      ? { label: '32', value: '32' }
      : matrixOptions[0];

  if (profile.value == 'CN')
    return options.includes('48')
      ? { label: '48', value: '48' }
      : matrixOptions[0];
}

export function defaultFilter2(filterOptions) {
  const options = filterOptions.map(({ value }) => value);
  return options.includes('NA')
    ? { label: 'NA', value: 'NA' }
    : filterOptions[0];
}

export function defaultSignatureSet(signatureSetOptions) {
  const options = signatureSetOptions.map(({ value }) => value);
  return options.includes('NA')
    ? { label: 'NA', value: 'NA' }
    : signatureSetOptions[0];
}

export function defaultSignatureSet2(signatureSetOptions) {
  const options = signatureSetOptions.map(({ value }) => value);
  return options.includes('COSMIC_v3.3_Signatures_GRCh38_SBS96')
    ? {
        label: 'COSMIC_v3.3_Signatures_GRCh38_SBS96',
        value: 'COSMIC_v3.3_Signatures_GRCh38_SBS96',
      }
    : signatureSetOptions[0];
}

export function defaultStrategy(strategyOptions) {
  const options = strategyOptions.map(({ value }) => value);
  return options.includes('NA')
    ? { label: 'NA', value: 'NA' }
    : strategyOptions[0];
}

export function defaultSignatureName(signatureNameOptions) {
  const options = signatureNameOptions.map(({ value }) => value);
  return options.includes('NA')
    ? { label: 'NA', value: 'NA' }
    : signatureNameOptions[0];
}

export function getJSON(path) {
  return fetch(`web/getFileS3`, {
    method: 'POST',
    headers: {
      Accept: 'application/json',
      'Content-Type': 'application/json',
    },
    body: JSON.stringify({
      path: path,
    }),
  })
    .then((res) => res.json())
    .then((data) => data);
}

export function getBlob(path) {
  return fetch(`web/getFileS3`, {
    method: 'POST',
    headers: {
      Accept: 'application/json',
      'Content-Type': 'application/json',
    },
    body: JSON.stringify({
      path: path,
    }),
  })
    .then((res) => res.blob())
    .then((data) => data);
}
