import { useRef, useEffect } from 'react';
import { renderToStaticMarkup } from 'react-dom/server';
import * as d3 from 'd3';
import cloneDeep from 'lodash/cloneDeep';
import { useRecoilValue } from 'recoil';
import { formState, graphDataSelector } from './treeLeaf.state';
import { groupBy } from './treeLeaf.utils';

export default function D3TreeLeaf({ id = 'treeleaf-plot', width = 1000, height = 1000, ...props }) {
  const plotRef = useRef(null);
  const form = useRecoilValue(formState);
  const graphData = useRecoilValue(graphDataSelector);
  const { hierarchy, attributes } = cloneDeep(graphData) || {};
  const plotData = {
    data: hierarchy,
    attributes: groupBy(attributes, 'Sample'),
    form: form,
  };
  const plotLayout = {
    id,
    width,
    height,
    radius: Math.min(width, height) / 2,
  };

  useEffect(() => {
    if (plotRef.current && plotData.data) {
      const plot = createForceDirectedTree(plotData, plotLayout);
      plotRef.current.replaceChildren(plot);
    }
  }, [plotRef, plotData, plotLayout]);

  return (
    <div className="border rounded p-3" {...props}>
      <div ref={plotRef} />
    </div>
  );
}

// Copyright 2022 Observable, Inc.
// Released under the ISC license.
// https://observablehq.com/@d3/force-directed-tree
function createForceDirectedTree(
  { data, attributes, form },
  {
    id,
    children, // if hierarchical data, given a d in data, returns its children
    tree = d3.tree, // layout algorithm (typically d3.tree or d3.cluster)
    separation = (a, b) => (a.parent == b.parent ? 1 : 1) / a.depth,
    sort, // how to sort nodes prior to layout (e.g., (a, b) => d3.descending(a.height, b.height))
    label = (d) => d.name, // given a node d, returns the display name
    title, // given a node d, returns its hover text
    width = 640, // outer width, in pixels
    height = 400, // outer height, in pixels
    margin = 0, // shorthand for margins
    marginTop = margin, // top margin, in pixels
    marginRight = margin, // right margin, in pixels
    marginBottom = margin, // bottom margin, in pixels
    marginLeft = margin, // left margin, in pixels
    radius = Math.min(
      width - marginLeft - marginRight,
      height - marginTop - marginBottom
    ) / 2, // outer radius
    r = d3.scaleSqrt().range([5, 20]), // radius of nodes
    padding = 1, // horizontal padding for first and last column
    fill = form.color.continuous
      ? d3.scaleSequential(d3.interpolateViridis)
      : d3.scaleOrdinal(d3.schemeCategory10), // fill for nodes
    fillOpacity, // fill opacity for nodes
    stroke = '#555', // stroke for links
    strokeWidth = 1.5, // stroke width for links
    strokeOpacity = 0.4, // stroke opacity for links
    strokeLinejoin, // stroke line join for links
    strokeLinecap, // stroke line cap for links
    halo = '#fff', // color of label halo
    haloWidth = 3, // padding around the labels
    scale = 2,
    highlightQuery = null,
  }
) {
  // gather range of attributes
  const mutations = Object.values(attributes).map((e) => e.Mutations);
  const mutationMin = d3.min(mutations);
  const mutationMax = d3.max(mutations);

  const colorValues = Object.values(attributes).map((e) => e[form.color.value]);
  const colorMin = d3.min(colorValues);
  const colorMax = d3.max(colorValues);
  const colorFill = form.color.continuous
    ? fill.domain([colorMin, colorMax])
    : fill.domain(colorValues);

  const treeData = d3.hierarchy(data, children);

  // Sort the nodes.
  if (sort != null) treeData.sort(sort);

  separation = (a, b) => (a.parent == b.parent ? 1 : 2) / a.depth;

  // Compute the layout.
  const root = d3
    .tree()
    .size([2 * Math.PI, radius])
    .separation(separation)(treeData);

  treeData.each((d) => {
    let angle = d.x;
    let distance = d.y;
    d.x = distance * Math.cos(angle) * scale;
    d.y = distance * Math.sin(angle) * scale;
  });

  // Compute labels and titles.
  const nodes = root.descendants();
  const links = root.links();

  const nodeRadius = ({ data, children }) =>
    !children && data.name && attributes[data.name]
      ? r.domain([mutationMin, mutationMax])(attributes[data.name].Mutations)
      : 0;

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
    .force(
      'charge',
      d3.forceManyBody().strength((d) => (d.children ? -15 : 15))
    )
    // .force("center", d3.forceCenter().strength(1))
    .force('x', d3.forceX().strength(0.005))
    .force('y', d3.forceY().strength(0.005))
    .force('collision', d3.forceCollide().radius(nodeRadius));

  // simulation.stop();
  // simulation.tick(10);

  const zoom = d3.zoom().on('zoom', zoomed);

  const container = d3.create('div');
  const svg = container
    .append('svg')
    .attr('id', id)
    .attr('viewBox', [-marginLeft - radius, -marginTop - radius, width, height])
    .attr('width', width)
    .attr('height', height)
    .attr(
      'style',
      'max-width: 100%; height: auto; height: intrinsic; width: 100%;'
    )
    .attr('font-family', 'sans-serif')
    .attr('font-size', 10)
    .call(zoom);

  // add lines
  const link = svg
    .append('g')
    .attr('id', 'treeleaf-links')
    .attr('fill', 'none')
    .attr('stroke', stroke)
    .attr('stroke-opacity', strokeOpacity)
    .attr('stroke-linecap', strokeLinecap)
    .attr('stroke-linejoin', strokeLinejoin)
    .attr('stroke-width', strokeWidth)
    .selectAll('path')
    .data(links)
    .join('line')
    .attr('x1', (d) => d.source.x)
    .attr('y1', (d) => d.source.y)
    .attr('x2', (d) => d.target.x)
    .attr('y2', (d) => d.target.y);

  const searchValues = form.searchSamples?.map((s) => s.value) || [];
  const highlightedColor = 'yellow';

  function getNodeColor({ data }) {
    if (data.name && searchValues.includes(data.name[0])) {
      return highlightedColor;
    }

    return data.name && attributes[data.name]
      ? colorFill(attributes[data.name][form.color.value])
      : stroke;
  }

  function getNodeOpacity({ data }) {
    if (
      !searchValues?.length ||
      (data.name && searchValues.includes(data.name[0]))
    ) {
      return 1;
    }

    return 0.8;
  }

  const node = svg
    .append('g')
    .attr('id', 'treeleaf-nodes')
    .attr('fill', '#fff')
    .attr('stroke', '#000')
    .attr('stroke-width', 1.5)
    .selectAll('circle')
    .data(nodes)
    .join('circle')
    .attr('id', (d) => d.data.name[0])
    .attr('fill', getNodeColor)
    .attr('opacity', getNodeOpacity)
    .attr('stroke', stroke)
    .attr('stroke-width', strokeWidth)
    .attr('r', nodeRadius)
    .attr('cx', (d) => d.x)
    .attr('cy', (d) => d.y)
    .on('mouseover', mouseover)
    .on('mousemove', mousemove)
    .on('mouseleave', mouseleave);

  simulation.on('tick', () => {
    link
      .attr('x1', (d) => d.source.x)
      .attr('y1', (d) => d.source.y)
      .attr('x2', (d) => d.target.x)
      .attr('y2', (d) => d.target.y);

    node.attr('cx', (d) => d.x).attr('cy', (d) => d.y);
  });

  if (title != null) node.append('title').text((d) => title(d.data, d));

  function zoomed({ transform }) {
    d3.selectAll('#treeleaf-links,#treeleaf-nodes').attr(
      'transform',
      transform
    );
  }

  // add tooltips
  const tooltip = container
    .append('div')
    .style('position', 'absolute')
    .style('visibility', 'hidden')
    .attr('class', 'bg-light border rounded p-1');

  function mouseover() {
    tooltip.style('visibility', 'visible');
  }

  function mouseleave() {
    tooltip.style('visibility', 'hidden');
  }

  function mousemove(e, d) {
    const sample = d.data.name;
    const data = attributes[sample];
    // console.log(data);
    const content = (
      <div className="text-start">
        <div>Sample: {sample}</div>
        <div>Cancer Type: {data.Cancer_Type}</div>
        <div>Cosine Similarity: {data.Cosine_similarity}</div>
        <div>Mutations: {data.Mutations}</div>
      </div>
    );
    tooltip
      .html(renderToStaticMarkup(content))
      // .html('Sample: ' + sample || 'none')
      .style('left', e.pageX - 30 + 'px')
      .style('top', e.pageY - 320 + 'px');
  }

  const legendParams = {
    svg,
    color: colorFill,
    title: form.color.label,
    x: width/2 - 40, 
    y: -height/2 + 20
  };

  form.color.continuous
    ? continuousLegend(legendParams)
    : categoricalLegend(legendParams);

  return container.node();
}

function continuousLegend({
  svg,
  color,
  title = '',
  x = 10,
  y = 10,
  width = 20,
  height = 200,
} = {}) {
  // define gradient
  const defs = svg.append('defs');
  const linearGradient = defs
    .append('linearGradient')
    .attr('id', 'linear-gradient')
    .attr('gradientTransform', 'rotate(90)');
  linearGradient
    .selectAll('stop')
    .data(
      color
        .ticks()
        .reverse()
        .map((t, i, n) => ({
          offset: `${(100 * i) / n.length}%`,
          color: color(t),
        }))
    )
    .enter()
    .append('stop')
    .attr('offset', (d) => d.offset)
    .attr('stop-color', (d) => d.color);

  // add legend rect
  const group = svg
    .append('g')
    .attr('id', 'treeleaf-legend')
    .attr('transform', `translate(${x}, ${y})`);

  const legendBar = group
    .append('g')
    .attr('transform', `translate(0, 15)`)

  legendBar
    .append('rect')
    .attr('width', width)
    .attr('height', height)
    .style('fill', 'url(#linear-gradient)');

  // add legend axis labels
  const axisScale = d3.scaleLinear().domain(color.domain()).range([height, 0]);
  const axisLeft = (g) => g.call(d3.axisLeft(axisScale));
  legendBar.call(axisLeft);

  // add title
  group
    .append('text')
    .attr('x', width)
    .attr('y', 0)
    .attr('fill', 'black')
    .attr('text-anchor', 'end')
    .attr('font-weight', 'bold')
    .attr('class', 'title')    
    .text(title);

  return group;
}

function categoricalLegend({
  svg,
  color,
  title,
  x = 10,
  y = 10,
  size = 15,
} = {}) {
  const keys = color.domain();

  // add legend group
  const group = svg
    .append('g')
    .attr('id', 'treeleaf-legend')
    .attr('transform', `translate(${x}, ${y})`);

  // add legend bar
  const legendBar = group
    .append('g')
    .attr('transform', `translate(0, 10)`)

  legendBar
    .selectAll('treeleaf-legend-colors')
    .data(keys)
    .enter()
    .append('rect')
    .attr('x', 0)
    .attr('y', (d, i) => i * size)
    .attr('width', size)
    .attr('height', size)
    .style('fill', color);

  legendBar
    .selectAll('treeleaf-legend-labels')
    .data(keys)
    .enter()
    .append('text')
    .attr('x', -5)
    .attr('y', (d, i) => i * (size) + size / 2)
    .attr('fill', 'black')
    .attr('text-anchor', 'end')
    .style('alignment-baseline', 'middle')
    .text((d) => d);

  // add title
  group
    .append('text')
    .attr('x', size)
    .attr('y', 0)
    .attr('fill', 'black')
    .attr('text-anchor', 'end')
    .attr('font-weight', 'bold')
    .attr('class', 'title')    
    .text(title);    

  return group;
}
