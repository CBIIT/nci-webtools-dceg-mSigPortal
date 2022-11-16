import { useRef, useEffect } from 'react';
import { renderToStaticMarkup } from 'react-dom/server';
import * as d3 from 'd3';
import * as htl from 'htl';
import cloneDeep from 'lodash/cloneDeep';
import { useRecoilValue } from 'recoil';
import { getGraphData, formState } from './treeLeaf.state';

export default function D3TreeLeaf({ width = 1000, height = 1000, ...props }) {
  const { hierarchy, attributes } = cloneDeep(useRecoilValue(getGraphData));
  const form = useRecoilValue(formState);

  const plotRef = useRef(null);
  const colorLegendRef = useRef(null);
  // const mutationsLegendRef = useRef(null);

  const groupAttributesBySample = attributes.reduce(
    (map, e) => ((map[e.Sample] = e), map),
    {}
  );

  useEffect(() => {
    if (plotRef.current && hierarchy) {
      console.log({hierarchy})
      const [plot, colorLegend, mutationsLegend] = createForceDirectedTree(
        hierarchy,
        groupAttributesBySample,
        form,
        {
          width,
          height,
          radius: Math.min(width, height) / 2,
        }
      );

      plotRef.current.replaceChildren(plot);
      colorLegendRef.current.replaceChildren(colorLegend);
      // mutationsLegendRef.current.replaceChildren(mutationsLegend);
    }
  }, [hierarchy, form, plotRef, width, height]);

  return (
    <div className="border rounded p-3" {...props}>
      <div ref={colorLegendRef} />
      {/* <div ref={mutationsLegendRef} /> */}
      <div ref={plotRef} />
    </div>
  );
}

// Copyright 2022 Observable, Inc.
// Released under the ISC license.
// https://observablehq.com/@d3/force-directed-tree
function createForceDirectedTree(
  data,
  attributes,
  form,
  {
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
  const root = d3.tree()
    .size([2 * Math.PI, radius])
    .separation(separation)(treeData);

  treeData.each(d => {
    let angle = d.x;
    let distance = d.y;
    d.x = distance * Math.cos(angle) * scale;
    d.y = distance * Math.sin(angle) * scale;
  });

  // Compute labels and titles.
  const nodes = root.descendants();
  const links = root.links();
  console.log(nodes, links);

  const nodeRadius = ({ data, children }) =>
  !children && data.name && attributes[data.name]
    ? r.domain([mutationMin, mutationMax])(attributes[data.name].Mutations)
    : 0;


  const simulation = d3
    .forceSimulation(nodes)
    .force(
      "link",
      d3
        .forceLink(links)
        .id((d) => d.id)
        .distance(1)
        .strength(1)
    )
    .force("charge", d3.forceManyBody().strength(d => { return d.children ? -15 : 15 }))
    // .force("center", d3.forceCenter().strength(1))
    .force("x", d3.forceX().strength(0.005))
    .force("y", d3.forceY().strength(0.005))
    .force('collision', d3.forceCollide().radius(nodeRadius))

  // simulation.stop();
  // simulation.tick(10);

  const zoom = d3
    .zoom()
    .on('zoom', zoomed)

  const container = d3.create('div');
  const svg = container
    .append('svg')
    .attr('id', 'plot')
    .attr('viewBox', [-marginLeft - radius, -marginTop - radius, width, height])
    .attr('width', width)
    .attr('height', height)
    .attr('style', 'max-width: 100%; height: auto; height: intrinsic;')
    .attr('font-family', 'sans-serif')
    .attr('font-size', 10)
    .call(zoom);

  // add lines
  const link = svg
    .append('g')
    .attr('id', 'links')
    .attr('fill', 'none')
    .attr('stroke', stroke)
    .attr('stroke-opacity', strokeOpacity)
    .attr('stroke-linecap', strokeLinecap)
    .attr('stroke-linejoin', strokeLinejoin)
    .attr('stroke-width', strokeWidth)
    .selectAll('path')
    .data(links)
    .join('line')
    .attr("x1", (d) => d.source.x)
    .attr("y1", (d) => d.source.y)
    .attr("x2", (d) => d.target.x)
    .attr("y2", (d) => d.target.y);
  
  const searchValues = form.searchSamples?.map(s => s.value) || [];
  const highlightedColor = 'yellow';

  function getNodeColor({ data }) {
    if (data.name && searchValues.includes(data.name[0])) {
      return highlightedColor;
    }

    return data.name && attributes[data.name]
      ? colorFill(attributes[data.name][form.color.value])
      : stroke
  }

  function getNodeOpacity({ data }) {
    if (!searchValues?.length || (data.name && searchValues.includes(data.name[0]))) {
      return 1;
    }

    return 0.8;
  }

  const node = svg
    .append("g")
    .attr('id', 'nodes')
    .attr("fill", "#fff")
    .attr("stroke", "#000")
    .attr("stroke-width", 1.5)
  .selectAll("circle")
  .data(nodes)
  .join("circle")
    .attr('id', d => d.data.name[0])
    .attr('fill', getNodeColor)
    .attr('opacity', getNodeOpacity)
    .attr('stroke', stroke)
    .attr('stroke-width', strokeWidth)
    .attr('r', nodeRadius)
    .attr("cx", (d) => d.x)
    .attr("cy", (d) => d.y)
    .on('mouseover', mouseover)
    .on('mousemove', mousemove)
    .on('mouseleave', mouseleave)

    simulation.on("tick", () => {
      link
          .attr("x1", d => d.source.x)
          .attr("y1", d => d.source.y)
          .attr("x2", d => d.target.x)
          .attr("y2", d => d.target.y);
  
      node
          .attr("cx", d => d.x)
          .attr("cy", d => d.y);
    });

  if (title != null) node.append('title').text((d) => title(d.data, d));

  function zoomed({ transform }) {
    d3.selectAll('#plot g').attr('transform', transform);
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


  const colorLegend = form.color.continuous
    ? ContinuousLegend(colorFill, {
        title: form.color.label,
      })
    : CategoricalLegned(colorFill);

  return [container.node(), colorLegend];
}


// Copyright 2022 Observable, Inc.
// Released under the ISC license.
// https://observablehq.com/@d3/radial-tree
function createRadialTree(
  data,
  attributes,
  form,
  {
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

  // Compute the layout.
  const root = tree()
    .size([2 * Math.PI, radius])
    .separation(separation)(treeData);

  // Compute labels and titles.
  const nodes = root.descendants();

  const links = root.links();

  const zoom = d3
    .zoom()
    .scaleExtent([1, 10])
    .on('zoom', zoomed)
    .on('end', rescale);

  function zoomed({ transform }) {
    d3.selectAll('#plot').attr('transform', transform);
  }

  function rescale({ transform }) {
    node
      .select('circle')
      .attr('r', ({ data }) =>
        data.name && attributes[data.name]
          ? r.domain([mutationMin, mutationMax])(
              attributes[data.name].Mutations
            ) / transform.k
          : null
      )
      .attr('stroke-width', strokeWidth / transform.k);

    node.select('text').attr('font-size', `${0.75 / transform.k}rem`);

    link.attr('stroke-width', strokeWidth / transform.k);
  }

  const container = d3.create('div');
  const svg = container
    .append('svg')
    .attr('viewBox', [-marginLeft - radius, -marginTop - radius, width, height])
    .attr('width', width)
    .attr('height', height)
    .attr('style', 'max-width: 100%; height: auto; height: intrinsic;')
    .attr('font-family', 'sans-serif')
    .attr('font-size', 10)
    // .call(zoom);

  // add lines
  const link = svg
    .append('g')
    .attr('id', 'plot')
    .attr('fill', 'none')
    .attr('stroke', stroke)
    .attr('stroke-opacity', strokeOpacity)
    .attr('stroke-linecap', strokeLinecap)
    .attr('stroke-linejoin', strokeLinejoin)
    .attr('stroke-width', strokeWidth)
    .selectAll('path')
    .data(links)
    .join('path')
    .attr(
      'd',
      d3
        .linkRadial()
        .angle((d) => d.x)
        .radius((d) => d.y)
    );

  // add sample nodes
  const node = svg
    .append('g')
    .attr('id', 'plot')
    .selectAll('a')
    .data(nodes)
    .join('a')
    // .attr('xlink:href', link == null ? null : (d) => link(d.data, d))
    // .attr('target', link == null ? null : linkTarget)
    .attr(
      'transform',
      (d) => `rotate(${(d.x * 180) / Math.PI - 90}) translate(${d.y},0)`
    );

  node
    .append('circle')
    .attr('fill', ({ data }) =>
      data.name && attributes[data.name]
        ? colorFill(attributes[data.name][form.color.value])
        : stroke
    )
    .attr('stroke', stroke)
    .attr('stroke-width', strokeWidth)
    .attr('r', ({ data }) =>
      data.name && attributes[data.name]
        ? r.domain([mutationMin, mutationMax])(attributes[data.name].Mutations)
        : null
    )
    .on('mouseover', mouseover)
    .on('mousemove', mousemove)
    .on('mouseleave', mouseleave);

  if (form.showLabels) {
    node
      .append('text')
      .attr('transform', (d) => `rotate(${d.x >= Math.PI ? 180 : 0})`)
      .attr('dy', '0.32em')
      .attr('x', (d) => (d.x < Math.PI === !d.children ? 6 : -6))
      .attr('text-anchor', (d) =>
        d.x < Math.PI === !d.children ? 'start' : 'end'
      )
      .attr('paint-order', 'stroke')
      // .attr('stroke', halo)
      .attr('stroke-width', haloWidth)
      .attr('font-size', '.75rem')
      .text((d, i) => d.data.name);
  }

  if (title != null) node.append('title').text((d) => title(d.data, d));

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
    console.log(data);
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

  const colorLegend = form.color.continuous
    ? ContinuousLegend(colorFill, {
        title: form.color.label,
      })
    : CategoricalLegned(colorFill);

  return [container.node(), colorLegend];
}

// Copyright 2021, Observable Inc.
// Released under the ISC license.
// https://observablehq.com/@d3/color-legend
function ContinuousLegend(
  color,
  {
    title,
    tickSize = 6,
    width = 320,
    height = 44 + tickSize,
    marginTop = 18,
    marginRight = 0,
    marginBottom = 16 + tickSize,
    marginLeft = 0,
    ticks = width / 64,
    tickFormat,
    tickValues,
  } = {}
) {
  function ramp(color, n = 256) {
    const canvas = document.createElement('canvas');
    canvas.width = n;
    canvas.height = 1;
    const context = canvas.getContext('2d');
    for (let i = 0; i < n; ++i) {
      context.fillStyle = color(i / (n - 1));
      context.fillRect(i, 0, 1, 1);
    }
    return canvas;
  }

  const svg = d3
    .create('svg')
    .attr('width', width)
    .attr('height', height)
    .attr('viewBox', [0, 0, width, height])
    .style('overflow', 'visible')
    .style('display', 'block');

  let tickAdjust = (g) =>
    g.selectAll('.tick line').attr('y1', marginTop + marginBottom - height);
  let x;

  // Sequential
  if (color.interpolator) {
    x = Object.assign(
      color
        .copy()
        .interpolator(d3.interpolateRound(marginLeft, width - marginRight)),
      {
        range() {
          return [marginLeft, width - marginRight];
        },
      }
    );

    svg
      .append('image')
      .attr('x', marginLeft)
      .attr('y', marginTop)
      .attr('width', width - marginLeft - marginRight)
      .attr('height', height - marginTop - marginBottom)
      .attr('preserveAspectRatio', 'none')
      .attr('xlink:href', ramp(color.interpolator()).toDataURL());

    // scaleSequentialQuantile doesn’t implement ticks or tickFormat.
    if (!x.ticks) {
      if (tickValues === undefined) {
        const n = Math.round(ticks + 1);
        tickValues = d3
          .range(n)
          .map((i) => d3.quantile(color.domain(), i / (n - 1)));
      }
      if (typeof tickFormat !== 'function') {
        tickFormat = d3.format(tickFormat === undefined ? ',f' : tickFormat);
      }
    }
  }

  svg
    .append('g')
    .attr('id', `legend-${title.replace(/\W/g, '')}`)
    .attr('transform', `translate(0,${height - marginBottom})`)
    .call(
      d3
        .axisBottom(x)
        .ticks(ticks, typeof tickFormat === 'string' ? tickFormat : undefined)
        .tickFormat(typeof tickFormat === 'function' ? tickFormat : undefined)
        .tickSize(tickSize)
        .tickValues(tickValues)
    )
    .call(tickAdjust)
    .call((g) => g.select('.domain').remove())
    .call((g) =>
      g
        .append('text')
        .attr('x', marginLeft)
        .attr('y', marginTop + marginBottom - height - 6)
        .attr('fill', 'currentColor')
        .attr('text-anchor', 'start')
        .attr('font-weight', 'bold')
        .attr('class', 'title')
        .text(title)
    );

  return svg.node();
}
function CategoricalLegned(
  color,
  {
    columns = null,
    format,
    unknown: formatUnknown,
    swatchSize = 15,
    swatchWidth = swatchSize,
    swatchHeight = swatchSize,
    marginLeft = 0,
  } = {}
) {
  const id = `-swatches-${Math.random().toString(16).slice(2)}`;
  const unknown = formatUnknown == null ? undefined : color.unknown();
  const unknowns =
    unknown == null || unknown === d3.scaleImplicit ? [] : [unknown];
  const domain = color.domain().concat(unknowns);
  if (format === undefined) format = (x) => (x === unknown ? formatUnknown : x);

  function entity(character) {
    return `&#${character.charCodeAt(0).toString()};`;
  }

  if (columns !== null)
    return htl.html`<div style="display: flex; align-items: center; margin-left: ${+marginLeft}px; min-height: 33px; font: 10px sans-serif;">
  <style>

.${id}-item {
  break-inside: avoid;
  display: flex;
  align-items: center;
  padding-bottom: 1px;
}

.${id}-label {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  max-width: calc(100% - ${+swatchWidth}px - 0.5em);
}

.${id}-swatch {
  width: ${+swatchWidth}px;
  height: ${+swatchHeight}px;
  margin: 0 0.5em 0 0;
}

  </style>
  <div style=${{ width: '100%', columns }}>${domain.map((value) => {
      const label = `${format(value)}`;
      return htl.html`<div class=${id}-item>
      <div class=${id}-swatch style=${{ background: color(value) }}></div>
      <div class=${id}-label title=${label}>${label}</div>
    </div>`;
    })}
  </div>
</div>`;

  return htl.html`<div style="display: flex; align-items: center; min-height: 33px; margin-left: ${+marginLeft}px; font: 10px sans-serif;">
  <style>

.${id} {
  display: inline-flex;
  align-items: center;
  margin-right: 1em;
}

.${id}::before {
  content: "";
  width: ${+swatchWidth}px;
  height: ${+swatchHeight}px;
  margin-right: 0.5em;
  background: var(--color);
}

  </style>
  <div>${domain.map(
    (value) =>
      htl.html`<span class="${id}" style="--color: ${color(value)}">${format(
        value
      )}</span>`
  )}</div>`;
}
