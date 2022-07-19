// organize seqmatrix data into a format suitable for graphing in plotly
const regexMap = {
  SBS96: /\[(.*)\]/,
  DBS78: /^(.{2})/,
  ID83: /^(.{7})/,
};

function baseSubstitution(data, profile, matrix) {
  const groupByMutation = data.reduce((acc, e, i) => {
    const mutationRegex = regexMap[profile + matrix];
    const mutation = e.mutationType.match(mutationRegex)[1];

    acc[mutation] = acc[mutation] ? [...acc[mutation], e] : [e];
    return acc;
  }, {});

  const transform = Object.entries(groupByMutation).map(([mutation, data]) => ({
    mutation,
    data,
  }));

  return transform;
}

function indel(data, matrix) {
  // return data;
  const groupByIndel = data.reduce((acc, e, i) => {
    const indel = e.mutationType.match(regexMap['ID' + matrix])[1];

    acc[indel] = acc[indel] ? [...acc[indel], e] : [e];
    return acc;
  }, {});

  const transform = Object.entries(groupByIndel).map(([indel, data]) => ({
    indel,
    data,
  }));

  return transform;
}

module.exports = { baseSubstitution, indel };
