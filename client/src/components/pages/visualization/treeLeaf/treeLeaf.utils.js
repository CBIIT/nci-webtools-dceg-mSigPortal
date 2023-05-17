import { quadtree } from 'd3-quadtree';

export function exportSvg(selector, filename) {
  const svgString = new XMLSerializer().serializeToString(
    document.querySelector(selector)
  );
  const svgBlob = new Blob([svgString], { type: 'image/svg+xml' });
  const svgUrl = URL.createObjectURL(svgBlob);
  const downloadLink = document.createElement('a');
  downloadLink.href = svgUrl;
  downloadLink.download = filename;
  document.body.appendChild(downloadLink);
  downloadLink.click();
  document.body.removeChild(downloadLink);
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

export function createIdGenerator() {
  let start = 0;
  return () => ++start;
}

export function assignParents(node, getId) {
  node.id = getId();
  let children = Array.isArray(node.children) ? node.children : [node.children];
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
  const nodes = getAllChildren(node); // || getLeafNodes(node);// || getAllChildren(node);
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

export function createPromiseWorker(url, options) {
  const worker = new Worker(url, options);
  return {
    submit: (params) =>
      new Promise((resolve, reject) => {
        worker.onmessage = (event) => {
          resolve(event.data);
        };
        worker.onerror = (error) => {
          reject(error);
        };
        worker.postMessage(params);
      }),
    terminate: () => worker.terminate(),
  };
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
        const distanceFactor = 0.003;
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

export function forceHierarchical(attractionMatrix) {
  let nodes,
    node,
    random,
    alpha,
    strength = () => -30,
    strengths,
    distanceMin2 = 1,
    distanceMax2 = Infinity,
    theta2 = 0.81,
    x = (d) => d.x,
    y = (d) => d.y;

  function jiggle(random) {
    return (random() - 0.5) * 1e-6;
  }

  function force(_) {
    let i,
      n = nodes.length,
      tree = quadtree(nodes, x, y).visitAfter(accumulate);
    for (alpha = _, i = 0; i < n; ++i) {
      node = nodes[i];
      tree.visit(apply);
    }
  }

  function initialize() {
    if (!nodes) return;
    let i,
      n = nodes.length,
      node;
    strengths = new Array(n);
    for (i = 0; i < n; ++i) {
      node = nodes[i];
      strengths[node.index] = +strength(node, i, nodes);
    }
  }

  function accumulate(quad) {
    let strength = 0,
      q,
      c,
      weight = 0,
      x,
      y,
      i;

    // For internal nodes, accumulate forces from child quadrants.
    if (quad.length) {
      for (x = y = i = 0; i < 4; ++i) {
        if ((q = quad[i]) && (c = Math.abs(q.value))) {
          strength += q.value;
          weight += c;
          x += c * q.x;
          y += c * q.y;
        }
      }
      quad.x = x / weight;
      quad.y = y / weight;
    }

    // For leaf nodes, accumulate forces from coincident quadrants.
    else {
      q = quad;
      q.x = q.data.x;
      q.y = q.data.y;
      do strength += strengths[q.data.index];
      while ((q = q.next));
    }

    quad.value = strength;
  }

  function apply(quad, x1, _, x2) {
    if (!quad.value) return true;

    let x = quad.x - node.x,
      y = quad.y - node.y,
      w = x2 - x1,
      l = x * x + y * y;

    // Apply the Barnes-Hut approximation if possible.
    // Limit forces for very close nodes; randomize direction if coincident.
    if ((w * w) / theta2 < l) {
      if (l < distanceMax2) {
        if (x === 0) {
          x = jiggle(random);
          l += x * x;
        }
        if (y === 0) {
          y = jiggle(random);
          l += y * y;
        }
        if (l < distanceMin2) l = Math.sqrt(distanceMin2 * l);
        node.vx += (x * quad.value * alpha) / l;
        node.vy += (y * quad.value * alpha) / l;
      }
      return true;
    }

    // Otherwise, process points directly.
    else if (quad.length || l >= distanceMax2) return;

    // Limit forces for very close nodes; randomize direction if coincident.
    if (quad.data !== node || quad.next) {
      if (x === 0) {
        x = jiggle(random);
        l += x * x;
      }
      if (y === 0) {
        y = jiggle(random);
        l += y * y;
      }
      if (l < distanceMin2) l = Math.sqrt(distanceMin2 * l);
    }

    do
      if (quad.data !== node) {
        w = (strengths[quad.data.index] * alpha) / l;
        node.vx += x * w;
        node.vy += y * w;
      }
    while ((quad = quad.next));
  }

  force.initialize = function (_nodes, _random) {
    nodes = _nodes;
    random = _random;
    initialize();
  };

  force.strength = function (_) {
    return arguments.length
      ? ((strength = typeof _ === 'function' ? _ : () => +_),
        initialize(),
        force)
      : strength;
  };

  force.distanceMin = function (_) {
    return arguments.length
      ? ((distanceMin2 = _ * _), force)
      : Math.sqrt(distanceMin2);
  };

  force.distanceMax = function (_) {
    return arguments.length
      ? ((distanceMax2 = _ * _), force)
      : Math.sqrt(distanceMax2);
  };

  force.theta = function (_) {
    return arguments.length ? ((theta2 = _ * _), force) : Math.sqrt(theta2);
  };

  return force;
}
