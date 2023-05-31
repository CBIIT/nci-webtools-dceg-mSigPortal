import * as d3 from 'https://cdn.jsdelivr.net/npm/d3@7/+esm';

onmessage = (event) => {
  postMessage(createForceDirectedTree(event.data));
};

function createForceDirectedTree({
  data,
  attributes,
  radius = 300, // outer radius
}) {
  // gather range of attributes
  const separation = (a, b) => (a.parent === b.parent ? 1 : 2) / a.depth;
  simplifyTree(data);
  assignParents(data, createIdGenerator());

  // const attributeLookup = groupBy(attributes, 'name')
  const mutations = Object.values(attributes).map((e) => e.Mutations);
  const mutationMin = d3.min(mutations);
  const mutationMax = d3.max(mutations);
  const treeData = d3.hierarchy(data);
  const scale = Math.max(0.5, Math.log10(Object.keys(attributes).length) / 10);

  // Compute the layout.
  const root = d3
    .tree()
    .size([2 * Math.PI, radius])
    .separation(separation)(treeData);

  const r = d3.scaleSqrt().range([2, 10]);
  treeData.each((d) => {
    let angle = d.x;
    let distance = d.y;
    d.x = distance * Math.cos(angle) * scale;
    d.y = distance * Math.sin(angle) * scale;
    d.r = !d.children && d.data.name && attributes[d.data.name] && attributes[d.data.name].Mutations
    ? r.domain([mutationMin, mutationMax])(attributes[d.data.name].Mutations)
    : 0
  });

  const distanceMatrix = getAttractionMatrix(data, (node) => node.id);

  // Compute labels and titles.
  const nodes = root.descendants();
  const links = root.links();

  const simulation = d3
    .forceSimulation(nodes)
    .force('link', d3.forceLink(links).id((d) => d.id).distance(1).strength(1))
    .force('hierarchical', hierarchicalForce(distanceMatrix, nodes))
    // .force('manyBody', d3.forceManyBody().strength(d => d.children ? -5 : -15))
    // .force('manyBody', d3.forceManyBody().strength(-1))
    .force('center', d3.forceCenter().strength(1))
    .force('x', d3.forceX().strength(0.005))
    .force('y', d3.forceY().strength(0.005))
    .force('collision', d3.forceCollide().radius(d => d.r * 1.2));
  
  simulation.stop();
  simulation.tick(120);

  return { nodes, links };
}

export function hierarchicalForce(distanceMatrix, nodes) {
  return function (alpha) {
    for (var i = 0; i < nodes.length; i++) {
      for (var j = i + 1; j < nodes.length; j++) {
        const a = nodes[i];
        const b = nodes[j];
        const idA = a?.data?.id;
        const idB = b?.data?.id;

        // only include leaf nodes in the comparison
        if (a.depth === 0 || b.depth === 0 || !a.data.name || !b.data.name)
          continue;
        const distance =
          distanceMatrix[idA]?.[idB] || distanceMatrix[idB]?.[idA];
        if (!distance) continue;

        // get the angle between a and b
        const dx = b.x - a.x;
        const dy = b.y - a.y;
        const angle = Math.atan2(dy, dx);

        // move each node away from each other, multiplied by the scaled distance
        const distanceFactor = 0.001;
        const delta =
          (distance - 4 - Math.floor(Math.log2(nodes.length))) *
          distanceFactor *
          alpha ** 2;
        a.x -= Math.cos(angle) * delta;
        a.y -= Math.sin(angle) * delta;
        b.x += Math.cos(angle) * delta;
        b.y += Math.sin(angle) * delta;
      }
    }
  };
}

export function createIdGenerator() {
  let start = 0;
  return () => ++start;
}

export function getChildren(node) {
  return Array.isArray(node.children) ? node.children : [node.children].filter(Boolean);
}

// ensures nodes with only one child and one parent are removed
export function simplifyTree(node) {
  let children = getChildren(node);
  // debugger;
  if (children.length === 1) {
    console.log('simplifying 1 child', node)
    let grandchildren = getChildren(children[0]);
    if (grandchildren.length === 1) {
      console.log('simplifying 1 grandchild', node)
      node.children = grandchildren;
    }
  }
  
  for (let child of getChildren(node)) {
    simplifyTree(child);
  }
}

export function assignParents(node, getId) {
  node.id = getId();
  let children = getChildren(node);
  for (let child of children) {
    if (child) {
      child.parent = node;
      assignParents(child, getId);
    }
  }
  return node;
}

export function getParents(node) {
  if (node.parents && node.parents.length) return node.parents;

  const parents = [];
  let parent = node.parent;
  while (parent) {
    parents.push(parent);
    parent = parent.parent;
  }
  node.parents = parents;
  return parents;
}

export function getLeafNodes(node) {
  // iterate through the tree and return all leaf nodes
  const leafNodes = [];
  const stack = [node];
  while (stack.length) {
    const node = stack.pop();
    if (node && node.children) {
      if (!Array.isArray(node.children)) {
        node.children = [node.children];
      }
      stack.push(...node.children);
    } else if (node) {
      leafNodes.push(node);
    }
  }
  return leafNodes;
}

export function getAllChildren(node, nodes = []) {
  nodes.push(node);
  if (node.children) {
    if (!Array.isArray(node.children)) {
      node.children = [node.children];
    }
    for (let child of node.children) {
      getAllChildren(child, nodes);
    }
  }
  return nodes;
}

/**
 * Returns the distance to the closest common parent of two nodes
 * @param {any} nodeA
 * @param {any} nodeB
 */
export function getClosestParentDistance(nodeA, nodeB) {
  const parentsA = getParents(nodeA);
  const parentsB = getParents(nodeB);
  let distanceA = 0;
  let distanceB = 0;

  // find the sum of the distances to the closest common parent
  for (let i = 0; i < parentsA.length; i++) {
    const parentA = parentsA[i];
    const indexB = parentsB.indexOf(parentA);
    if (indexB !== -1) {
      distanceA = i;
      distanceB = indexB;
      break;
    }
  }

  return distanceA + distanceB;
}

export function getAttractionMatrix(node, id = (node) => node.id) {
  let matrix = {};
  const nodes = getAllChildren(node);
  for (let i = 0; i < nodes.length; i++) {
    for (let j = i + 1; j < nodes.length; j++) {
      let a = nodes[i];
      let b = nodes[j];

      let distance = getClosestParentDistance(a, b);
      let idA = id(nodes[i]);
      let idB = id(nodes[j]);

      matrix[idA] = matrix[idA] || {};
      matrix[idB] = matrix[idB] || {};
      matrix[idA][idB] = distance;
      matrix[idB][idA] = distance;
    }
  }

  return matrix;
}


export function groupBy(array, key) {
  return array.reduce(
    (result, currentValue) => ({
      ...result,
      [currentValue[key]]: currentValue,
    }),
    {}
  );
}