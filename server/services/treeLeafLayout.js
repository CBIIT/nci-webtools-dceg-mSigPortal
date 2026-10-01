import * as d3 from 'd3';

export function createForceDirectedTree({
  data,
  mutations: mutationsBySample,
  radius = 300, // outer radius
}) {
  const separation = (a, b) => (a.parent === b.parent ? 1 : 2) / a.depth;

  const mutations = [...mutationsBySample.values()];
  const mutationMin = d3.min(mutations);
  const mutationMax = d3.max(mutations);
  const treeData = d3.hierarchy(data);
  const scale = Math.max(0.5, Math.log10(mutationsBySample.size) / 2);

  const root = d3
    .tree()
    .size([2 * Math.PI, radius])
    .separation(separation)(treeData);

  const r = d3
    .scaleLog()
    .domain([mutationMin, mutationMax].map((e) => Math.max(1, e)))
    .range([2, 10]);
  treeData.each((d) => {
    const angle = d.x;
    const distance = d.y;
    const mutations =
      !d.children && d.data.name && mutationsBySample.get(String(d.data.name));
    d.x = distance * Math.cos(angle) * scale;
    d.y = distance * Math.sin(angle) * scale;
    d.r = mutations ? r(mutations) : 0;
  });

  const nodes = root.descendants();
  const links = root.links();

  const simulation = d3
    .forceSimulation(nodes)
    .force(
      'link',
      d3
        .forceLink(links)
        .id((d) => d.id)
        .distance(1)
        .strength(1)
    )
    .force('hierarchical', hierarchicalForce(root, nodes))
    .force('center', d3.forceCenter().strength(1))
    .force('x', d3.forceX().strength(0.005))
    .force('y', d3.forceY().strength(0.005))
    .force(
      'collision',
      d3
        .forceCollide()
        .radius((d) => (d.children ? 1 : d.r * (d.r < 4 ? 1.5 : 1.25)))
    );

  simulation.stop();
  simulation.tick(120);

  // plain records so the d3 hierarchy (with parent/child cycles) isn't retained;
  // links reference nodes by index (d.index is set by forceSimulation)
  return {
    nodes: nodes.map((d) => ({
      x: d.x,
      y: d.y,
      r: d.r,
      name: String(d.data.name ?? ''),
    })),
    links: links.map(({ source, target }) => [source.index, target.index]),
  };
}

function hierarchicalForce(root, nodes) {
  // names are arrays (e.g. [""]), so internal nodes are included in the filter
  const pairNodes = nodes.filter((d) => d.depth !== 0 && d.data.name);
  const getLcaDepth = createLcaDepthQuery(root);
  const parentDepth = new Int32Array(pairNodes.length);
  const parentPosition = new Int32Array(pairNodes.length);
  for (let i = 0; i < pairNodes.length; i++) {
    parentDepth[i] = pairNodes[i].parent.depth;
    parentPosition[i] = getLcaDepth.position(pairNodes[i].parent);
  }
  const distanceFactor = 0.001;
  const offset = 4 + Math.floor(Math.log2(nodes.length));

  return function (alpha) {
    const alpha2 = alpha ** 2;
    for (let i = 0; i < pairNodes.length; i++) {
      const a = pairNodes[i];
      for (let j = i + 1; j < pairNodes.length; j++) {
        // edges from each node's parent up to their closest shared ancestor
        const distance =
          parentDepth[i] +
          parentDepth[j] -
          2 * getLcaDepth(parentPosition[i], parentPosition[j]);
        if (!distance) continue;

        const b = pairNodes[j];
        const dx = b.x - a.x;
        const dy = b.y - a.y;
        const angle = Math.atan2(dy, dx);

        const delta = (distance - offset) * distanceFactor * alpha2;
        const moveX = Math.cos(angle) * delta;
        const moveY = Math.sin(angle) * delta;
        a.x -= moveX;
        a.y -= moveY;
        b.x += moveX;
        b.y += moveY;
      }
    }
  };
}

/**
 * Euler tour + sparse table: O(1) depth of the lowest common ancestor of two nodes,
 * addressed by their first position in the tour.
 */
function createLcaDepthQuery(root) {
  const firstPosition = new Map();
  const tour = [];
  const stack = [{ node: root, next: 0 }];
  firstPosition.set(root, 0);
  tour.push(root.depth);
  while (stack.length) {
    const top = stack[stack.length - 1];
    const children = top.node.children;
    if (children && top.next < children.length) {
      const child = children[top.next++];
      firstPosition.set(child, tour.length);
      tour.push(child.depth);
      stack.push({ node: child, next: 0 });
    } else {
      stack.pop();
      if (stack.length) tour.push(stack[stack.length - 1].node.depth);
    }
  }

  const size = tour.length;
  const log = new Uint8Array(size + 1);
  for (let i = 2; i <= size; i++) log[i] = log[i >> 1] + 1;
  const levels = log[size] + 1;
  const table = new Int32Array(levels * size);
  table.set(tour);
  for (let k = 1; k < levels; k++) {
    const span = 1 << (k - 1);
    const row = k * size;
    const prev = (k - 1) * size;
    for (let i = 0; i + (1 << k) <= size; i++) {
      table[row + i] = Math.min(table[prev + i], table[prev + i + span]);
    }
  }

  function query(a, b) {
    const left = a < b ? a : b;
    const right = a < b ? b : a;
    const k = log[right - left + 1];
    const row = k * size;
    return Math.min(table[row + left], table[row + right - (1 << k) + 1]);
  }
  query.position = (node) => firstPosition.get(node);
  return query;
}
