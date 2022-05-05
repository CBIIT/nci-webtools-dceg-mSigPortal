import { useRef, useEffect } from 'react';
import * as d3 from 'd3';
import cloneDeep from 'lodash/cloneDeep';
import { useRecoilValue } from 'recoil';
import { getGraphData } from './treeLeaf.state';

export default function D3TreeLeaf({ width = 1000, height = 1000, ...props }) {
  const data = cloneDeep(useRecoilValue(getGraphData));
  const nodeRef = useRef(null);

  useEffect(() => {
    if (nodeRef.current && data) {
      // const treeLeaf = createTreeHierarchy(cloneDeep(data), width, height);
      // const treeLeaf = createTreeGraph(cloneDeep(data), width, height);
      const treeLeaf = createRadialTree(
        data.hierarchy,
        data.attributes.reduce((map, e) => ((map[e.Sample] = e), map), {}),
        {
          label: (d) => '', //d.name,
          width,
          height,
        }
      );

      nodeRef.current.replaceChildren(treeLeaf);
    }
  }, [data, nodeRef, width, height]);

  return <div ref={nodeRef} {...props} />;
}

export function createTreeHierarchy(data, width, height) {
  const root = d3.hierarchy(data.hierarchy);
  const links = root.links();
  const nodes = root.descendants();

  console.log(nodes[0]);
  console.log(links[0]);

  const simulation = d3
    .forceSimulation(nodes)
    .force(
      'link',
      d3
        .forceLink(links)
        .id((d) => d.id)
        .distance(0)
        .strength(1)
    )
    .force(
      'charge',
      d3.forceManyBody().strength((d) => (d.children ? -1000 : -50))
    )
    .force('x', d3.forceX())
    .force('y', d3.forceY())
    .force(
      'collide',
      d3.forceCollide().radius((d) => {
        // if parent has many children, reduce collide strength
        if (d.parent) {
          if (d.parent.children.length > 10) {
            return 7;
          }
        }

        // increase collide strength
        if (d.children && d.depth > 2) {
          return 18;
        }
        return 18 * 2;
      })
    );

  const svg = d3
    .create('svg')
    .attr('viewBox', [-width / 2, -height / 2, width, height]);

  const link = svg
    .append('g')
    .attr('stroke', '#999')
    .attr('stroke-opacity', 0.6)
    .selectAll('line')
    .data(links)
    .join('line');

  const node = svg
    .append('g')
    .attr('fill', '#fff')
    .attr('stroke', '#000')
    .attr('stroke-width', 1.5)
    .selectAll('circle')
    .data(nodes)
    .join('circle')
    .attr('fill', (d) => (d.children ? null : '#000'))
    .attr('stroke', (d) => (d.children ? null : '#fff'))
    .attr('r', (d) => (d.children ? null : 10))
    .call(drag(simulation));

  simulation.on('tick', () => {
    link
      .attr('x1', (d) => d.source.x)
      .attr('y1', (d) => d.source.y)
      .attr('x2', (d) => d.target.x)
      .attr('y2', (d) => d.target.y);

    node.attr('cx', (d) => d.x).attr('cy', (d) => d.y);
  });

  function drag(simulation) {
    function dragstarted(event, d) {
      if (!event.active) simulation.alphaTarget(0.3).restart();
      d.fx = d.x;
      d.fy = d.y;
    }

    function dragged(event, d) {
      d.fx = event.x;
      d.fy = event.y;
    }

    function dragended(event, d) {
      if (!event.active) simulation.alphaTarget(0);
      d.fx = null;
      d.fy = null;
    }

    return d3
      .drag()
      .on('start', dragstarted)
      .on('drag', dragged)
      .on('end', dragended);
  }

  return svg.node();
}

export function createTreeGraph(data, width, height) {
  const { nodes, links } = data.graph;

  console.log(nodes[0]);
  console.log(links[0]);

  const simulation = d3
    .forceSimulation(nodes)
    .force('link', d3.forceLink(links).distance(0).strength(1))
    .force(
      'charge',
      d3
        .forceManyBody()
        .strength((d, i) => (d.name.includes('Node') ? -1000 : -200))
    )
    .force('x', d3.forceX())
    .force('y', d3.forceY());

  const svg = d3
    .create('svg')
    .attr('viewBox', [-width / 2, -height / 2, width, height]);

  const link = svg
    .append('g')
    .attr('stroke', '#999')
    .attr('stroke-opacity', 0.6)
    .selectAll('line')
    .data(links)
    .join('line');

  const node = svg
    .append('g')
    .attr('fill', '#fff')
    .attr('stroke', '#000')
    .attr('stroke-width', 1.5)
    .selectAll('circle')
    .data(nodes)
    .join('circle')
    .attr('fill', (d) => (d.name.includes('Node') ? null : '#000'))
    .attr('stroke', (d) => (d.name.includes('Node') ? null : '#fff'))
    .attr('r', (d) => (d.name.includes('Node') ? null : 10))
    .call(drag(simulation));

  simulation.on('tick', () => {
    link
      .attr('x1', (d) => d.source.x)
      .attr('y1', (d) => d.source.y)
      .attr('x2', (d) => d.target.x)
      .attr('y2', (d) => d.target.y);

    node.attr('cx', (d) => d.x).attr('cy', (d) => d.y);
  });

  function drag(simulation) {
    function dragstarted(event, d) {
      if (!event.active) simulation.alphaTarget(0.3).restart();
      d.fx = d.x;
      d.fy = d.y;
    }

    function dragged(event, d) {
      d.fx = event.x;
      d.fy = event.y;
    }

    function dragended(event, d) {
      if (!event.active) simulation.alphaTarget(0);
      d.fx = null;
      d.fy = null;
    }

    return d3
      .drag()
      .on('start', dragstarted)
      .on('drag', dragged)
      .on('end', dragended);
  }

  return svg.node();
}

// Copyright 2022 Observable, Inc.
// Released under the ISC license.
// https://observablehq.com/@d3/radial-tree
function createRadialTree(
  data,
  attributes,
  {
    // data is either tabular (array of objects) or hierarchy (nested objects)
    path, // as an alternative to id and parentId, returns an array identifier, imputing internal nodes
    id = Array.isArray(data) ? (d) => d.id : null, // if tabular data, given a d in data, returns a unique identifier (string)
    parentId = Array.isArray(data) ? (d) => d.parentId : null, // if tabular data, given a node d, returns its parent’s identifier
    children, // if hierarchical data, given a d in data, returns its children
    tree = d3.tree, // layout algorithm (typically d3.tree or d3.cluster)
    separation = tree === d3.tree
      ? (a, b) => (a.parent == b.parent ? 1 : 2) / a.depth
      : (a, b) => (a.parent == b.parent ? 1 : 2),
    sort, // how to sort nodes prior to layout (e.g., (a, b) => d3.descending(a.height, b.height))
    label, // given a node d, returns the display name
    title, // given a node d, returns its hover text
    link, // given a node d, its link (if any)
    linkTarget = '_blank', // the target attribute for links (if any)
    width = 640, // outer width, in pixels
    height = 400, // outer height, in pixels
    margin = 60, // shorthand for margins
    marginTop = margin, // top margin, in pixels
    marginRight = margin, // right margin, in pixels
    marginBottom = margin, // bottom margin, in pixels
    marginLeft = margin, // left margin, in pixels
    radius = Math.min(
      width - marginLeft - marginRight,
      height - marginTop - marginBottom
    ) / 2, // outer radius
    r = d3
      .scaleLinear()
      .domain([
        0,
        d3.max(Object.values(attributes).map(({ Mutations }) => Mutations)),
      ])
      .range([4, 20]), // radius of nodes
    padding = 1, // horizontal padding for first and last column
    fill = d3.scaleSequential(d3.interpolateRdYlGn), // fill for nodes
    fillOpacity, // fill opacity for nodes
    stroke = '#555', // stroke for links
    strokeWidth = 1.5, // stroke width for links
    strokeOpacity = 0.4, // stroke opacity for links
    strokeLinejoin, // stroke line join for links
    strokeLinecap, // stroke line cap for links
    halo = '#fff', // color of label halo
    haloWidth = 3, // padding around the labels
  }
) {
  console.log(data);
  console.log(attributes['SP117655']);
  console.log(r(attributes['SP117655'].Mutations));
  console.log(typeof fill(attributes['SP117655'].Cosine_similarity));
  // If id and parentId options are specified, or the path option, use d3.stratify
  // to convert tabular data to a hierarchy; otherwise we assume that the data is
  // specified as an object {children} with nested objects (a.k.a. the “flare.json”
  // format), and use d3.hierarchy.
  const root =
    path != null
      ? d3.stratify().path(path)(data)
      : id != null || parentId != null
      ? d3.stratify().id(id).parentId(parentId)(data)
      : d3.hierarchy(data, children);

  // Sort the nodes.
  if (sort != null) root.sort(sort);

  // Compute labels and titles.
  const descendants = root.descendants();
  const L = label == null ? null : descendants.map((d) => label(d.data, d));

  // Compute the layout.
  tree()
    .size([2 * Math.PI, radius])
    .separation(separation)(root);

  const svg = d3
    .create('svg')
    .attr('viewBox', [-marginLeft - radius, -marginTop - radius, width, height])
    .attr('width', width)
    .attr('height', height)
    .attr('style', 'max-width: 100%; height: auto; height: intrinsic;')
    .attr('font-family', 'sans-serif')
    .attr('font-size', 10);

  svg
    .append('g')
    .attr('fill', 'none')
    .attr('stroke', stroke)
    .attr('stroke-opacity', strokeOpacity)
    .attr('stroke-linecap', strokeLinecap)
    .attr('stroke-linejoin', strokeLinejoin)
    .attr('stroke-width', strokeWidth)
    .selectAll('path')
    .data(root.links())
    .join('path')
    .attr(
      'd',
      d3
        .linkRadial()
        .angle((d) => d.x)
        .radius((d) => d.y)
    );

  const node = svg
    .append('g')
    .selectAll('a')
    .data(root.descendants())
    .join('a')
    .attr('xlink:href', link == null ? null : (d) => link(d.data, d))
    .attr('target', link == null ? null : linkTarget)
    .attr(
      'transform',
      (d) => `rotate(${(d.x * 180) / Math.PI - 90}) translate(${d.y},0)`
    );

  // set size and color
  node
    .append('circle')
    .attr('fill', ({ data }) =>
      data.name && attributes[data.name]
        ? fill(attributes[data.name].Cosine_similarity)
        : 'stroke'
    )
    .attr('r', ({ data }) =>
      data.name && attributes[data.name]
        ? r(attributes[data.name].Mutations)
        : null
    );

  // add tooltips
  // node
  //   .selectAll('circle')
  //   .on('mouseover.tooltip', (e) => {
  //     const div = d3.select(e.target).append('div');

  //     div.attr('id', 'details');
  //     div.attr('class', 'tooltip');

  //     let rows = div
  //       .append('table')
  //       .selectAll('tr')
  //       .data(Object.keys(e))
  //       .enter()
  //       .append('tr');

  //     console.log(div);
  //     console.log(e);

  //     rows.append('th').text((key) => key);
  //     rows.append('td').text((key) => e[key]);
  //   })
  // .on('mousemove.tooltip', (e) => {
  //   let div = d3.select('div#details');

  //   // get height of tooltip
  //   let bbox = div.node().getBoundingClientRect();

  //   div.style('left', e.clientX + 'px');
  //   div.style('top', e.clientY - bbox.height + 'px');
  // })
  // .on('mouseleave.tooltip', (d) => {
  //   console.log('left');
  //   d3.selectAll('div#details').remove();
  // });

  if (title != null) node.append('title').text((d) => title(d.data, d));

  if (L)
    node
      .append('text')
      .attr('transform', (d) => `rotate(${d.x >= Math.PI ? 180 : 0})`)
      .attr('dy', '0.32em')
      .attr('x', (d) => (d.x < Math.PI === !d.children ? 6 : -6))
      .attr('text-anchor', (d) =>
        d.x < Math.PI === !d.children ? 'start' : 'end'
      )
      .attr('paint-order', 'stroke')
      .attr('stroke', halo)
      .attr('stroke-width', haloWidth)
      .text((d, i) => L[i]);

  return svg.node();
}
