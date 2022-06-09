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
