function SBS(data) {
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

function DBS(data) {
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

function ID(data) {
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
