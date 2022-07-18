// organize seqmatrix data into a format suitable for graphing in plotly
const regexMap = {
  SBS96: /\[(.*)\]/,
  DBS78: /^(.{2})/,
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
  const groupByMutation = data.reduce((acc, e, i) => {
    const mutationRegex = /\[(.*)\]/;
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

module.exports = { baseSubstitution, indel };
