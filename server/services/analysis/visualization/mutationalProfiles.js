function SBS(data, matrix) {
  const groupByBase = data.reduce((acc, e, i) => {
    const baseRegex = /\[(.*)\]/;
    const base = e.mutationType.match(baseRegex)[1];

    acc[base] = acc[base] ? [...acc[base], e] : [e];
    return acc;
  }, {});

  const transform = Object.entries(groupByBase).map(
    ([base, mutationTypes]) => ({
      base,
      mutationTypes,
    })
  );

  return transform;
}

function DBS(data, matrix) {
  const groupByBase = data.reduce((acc, e, i) => {
    const baseRegex = /^(.{2})/;
    const base = e.mutationType.match(baseRegex)[1];

    acc[base] = acc[base] ? [...acc[base], e] : [e];
    return acc;
  }, {});

  const transform = Object.entries(groupByBase).map(
    ([base, mutationTypes]) => ({
      base,
      mutationTypes,
    })
  );
  return transform;
}

function ID(data, matrix) {
  const groupByBase = data.reduce((acc, e, i) => {
    const baseRegex = /\[(.*)\]/;
    const base = e.mutationType.match(baseRegex)[1];

    acc[base] = acc[base] ? [...acc[base], e] : [e];
    return acc;
  }, {});

  const transform = Object.entries(groupByBase).map(
    ([base, mutationTypes]) => ({
      base,
      mutationTypes,
    })
  );
  return transform;
}

module.exports = { SBS, DBS, ID };
