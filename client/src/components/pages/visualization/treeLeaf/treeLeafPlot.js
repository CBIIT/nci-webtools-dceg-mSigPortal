import { useRef, useEffect } from 'react';
import { useSelector } from 'react-redux';
import * as d3 from 'd3';
import cloneDeep from 'lodash/cloneDeep';
import { useRecoilValue } from 'recoil';
import { formState, graphDataSelector } from './treeLeaf.state';
import {
  groupBy,
  forceHierarchical,
  getAttractionMatrix,
} from './treeLeaf.utils';

export default function D3TreeLeaf({
  id = 'treeleaf-plot',
  width = 1000,
  height = 1000,
  onSelect,
  ...props
}) {
  const plotRef = useRef(null);
  const form = useRecoilValue(formState);
  const publicForm = useSelector((state) => state.visualization.publicForm);
  const study = publicForm?.study?.value ?? 'PCAWG';
  const strategy = publicForm?.strategy?.value ?? 'WGS';
  const signatureSetName = 'COSMIC_v3_Signatures_GRCh37_SBS96';
  const profileMatrix = ['SBS96', 'DBS78', 'ID83'];
  const params = { study, strategy,  signatureSetName, profileMatrix, cancer: form?.cancerType?.value };
  const graphData = useRecoilValue(graphDataSelector(params));
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
    plotTitle: `${publicForm?.study?.label} - ${form.color.label}`,
  };
  const plotEvents = {
    onClick: onSelect,
  };


  useEffect(() => {
    if (plotRef.current && plotData.data) {
      const plot = createForceDirectedTree(plotData, plotLayout, plotEvents);
      plotRef.current.replaceChildren(plot);
    }
  }, [plotRef, plotData, plotLayout]);

  return (
    <div className="border rounded p-3 position-relative" {...props}>
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
    separation = (a, b) => (a.parent == b.parent ? 1 : 2) / a.depth,
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
      ? d3.scaleSequential(d3.interpolateBlues)
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
    plotTitle,
  },
  { onClick }
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

  const attractionMatrix = getAttractionMatrix();

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
      forceHierarchical(attractionMatrix).strength((d) =>
        d.children ? -15 : 15
      )
    )
    // .force("center", d3.forceCenter().strength(1))
    .force('x', d3.forceX().strength(0.005))
    .force('y', d3.forceY().strength(0.005))
    .force('collision', d3.forceCollide().radius(nodeRadius));

  // simulation.stop();
  // simulation.tick(40);

  const zoom = d3.zoom().on('zoom', zoomed);

  const viewBoxScale = 1.3;
  const container = d3.create('div');
  const svg = container
    .append('svg')
    .attr('id', id)
    .attr(
      'viewBox',
      [-marginLeft - radius, -marginTop - radius, width, height].map(
        (v) => v * viewBoxScale
      )
    )
    .attr('width', width)
    .attr('height', height)
    .attr(
      'style',
      'max-width: 100%; height: auto; height: intrinsic; width: 100%;'
    )
    .attr('font-family', 'sans-serif')
    .attr('font-size', 10)
    .call(zoom);

  // add tree container
  const treeGroup = svg.append('g').attr('transform', `translate(0, 20)`);

  // add lines
  const link = treeGroup
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
    if (data.name && searchValues.includes(data.name)) {
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

  const node = treeGroup
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
    .attr('class', 'c-pointer')
    .on('mouseover', mouseover)
    .on('mousemove', mousemove)
    .on('mouseleave', mouseleave)
    .on('click', click);

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
    tooltip
      .html(
        `<div class="text-start">
          <div>Sample: ${sample ?? 'Unavailable'}</div>
          <div>Cancer Type: ${data.Cancer_Type ?? 'Unavailable'}</div>
          <div>Cosine Similarity: ${
            data.Cosine_similarity ?? 'Unavailable'
          }</div>
          <div>Mutations: ${data.Mutations ?? 'Unavailable'}</div>
        </div>`
      )
      .style('left', e.layerX + 'px')
      .style('top', e.layerY + 'px');
  }

  function click(e, d) {
    const sample = d.data.name;
    const data = attributes[sample];
    if (typeof onClick === 'function') {
      onClick(data);
    }
  }

  // append title
  svg
    .append('g')
    .attr('id', 'treeleaf-title')
    .attr('transform', `translate(0, ${(-height / 2 + 10) * viewBoxScale})`)
    .append('text')
    .attr('x', 0)
    .attr('y', 20)
    .attr('fill', 'black')
    .attr('text-anchor', 'center')
    .style('font-weight', 'bold')
    .style('font-size', '24px')
    //.attr('class', 'h2')
    .text(plotTitle);

  // append legend
  const legendParams = {
    svg,
    color: colorFill,
    title: form.color.label,
    x: (width / 2 - 40) * viewBoxScale,
    y: (-height / 2 + 20) * viewBoxScale,
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
  width = 50,
  height = 300,
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

  const legendBar = group.append('g').attr('transform', `translate(0, 15)`);

  legendBar
    .append('rect')
    .attr('width', width)
    .attr('height', height)
    .style('fill', 'url(#linear-gradient)');

  // add legend axis labels
  const axisScale = d3.scaleLinear().domain(color.domain()).range([height, 0]);
  const axisLeft = (g) =>
    g.call(d3.axisLeft(axisScale)).style('font-size', '17px');
  legendBar.call(axisLeft);

  // add title
  group
    .append('text')
    .attr('x', width)
    .attr('y', 0)
    .attr('fill', 'black')
    .attr('text-anchor', 'end')
    .attr('font-weight', 'bold')
    .attr('class', 'h5 title')
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
  const legendBar = group.append('g').attr('transform', `translate(0, 10)`);

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
    .attr('y', (d, i) => i * size + size / 2)
    .attr('fill', 'black')
    .attr('text-anchor', 'end')
    .style('alignment-baseline', 'middle')
    .style('font-size', '17px')
    .text((d) => d ?? 'Unavailable');

  // add title
  group
    .append('text')
    .attr('x', size)
    .attr('y', 0)
    .attr('fill', 'black')
    .attr('text-anchor', 'end')
    .attr('font-weight', 'bold')
    .attr('class', 'h5 title')
    .text(title);

  return group;
}
