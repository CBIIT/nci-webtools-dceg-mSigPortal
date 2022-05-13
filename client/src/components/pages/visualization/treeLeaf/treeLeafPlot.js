import { useRef, useEffect } from 'react';
import { renderToStaticMarkup } from 'react-dom/server';
import * as d3 from 'd3';
import cloneDeep from 'lodash/cloneDeep';
import { useRecoilValue } from 'recoil';
import { getGraphData, formState } from './treeLeaf.state';

export default function D3TreeLeaf({ width = 1000, height = 1000, ...props }) {
  const { hierarchy, attributes } = cloneDeep(useRecoilValue(getGraphData));
  const form = useRecoilValue(formState);

  const plotRef = useRef(null);
  const csLegendRef = useRef(null);
  // const mutationsLegendRef = useRef(null);

  const groupAttributesBySample = attributes.reduce(
    (map, e) => ((map[e.Sample] = e), map),
    {}
  );

  useEffect(() => {
    if (plotRef.current && hierarchy) {
      const [plot, csLegend, mutationsLegend] = createRadialTree(
        hierarchy,
        groupAttributesBySample,
        {
          label: (d) => (form.showLabels ? d.name : ''),
          width,
          height,
        }
      );

      plotRef.current.replaceChildren(plot);
      csLegendRef.current.replaceChildren(csLegend);
      // mutationsLegendRef.current.replaceChildren(mutationsLegend);
    }
  }, [hierarchy, plotRef, width, height]);

  return (
    <div className="border rounded p-3" {...props}>
      <div ref={csLegendRef} />
      {/* <div ref={mutationsLegendRef} /> */}
      <div ref={plotRef} />
    </div>
  );
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
    separation = (a, b) => (a.parent == b.parent ? 1 : 2) / a.depth,
    sort, // how to sort nodes prior to layout (e.g., (a, b) => d3.descending(a.height, b.height))
    label, // given a node d, returns the display name
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
    r = d3.scaleSqrt().range([4, 20]), // radius of nodes
    padding = 1, // horizontal padding for first and last column
    fill = d3.scaleSequential(d3.interpolatePRGn), // fill for nodes
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
  const nodes = root.descendants();
  const L = label == null ? null : nodes.map((d) => label(d.data, d));

  const links = root.links();

  // Compute the layout.
  tree()
    .size([2 * Math.PI, radius])
    .separation(separation)(root);

  // const simulation = d3
  //   .forceSimulation(nodes)
  //   .force(
  //     'link',
  //     d3
  //       .forceLink(links)
  //       .id((d) => d.id)
  //       .distance(0)
  //       .strength(1)
  //   )
  //   .force('charge', d3.forceManyBody().strength(-10));

  const drag = (simulation) => {
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
  };

  const zoom = d3.zoom().scaleExtent([1, 5]).on('zoom', zoomed);

  function zoomed({ transform }) {
    d3.selectAll('#plot').attr('transform', transform);
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
    .call(zoom);

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

  // gather range of attributes
  const mutations = Object.values(attributes).map((e) => e.Mutations);
  const cs = Object.values(attributes).map((e) => e.Cosine_similarity);
  const mutationMin = d3.min(mutations);
  const mutationMax = d3.max(mutations);
  const csMin = d3.min(cs);
  const csMax = d3.max(cs);

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
    )
    .append('circle')
    .attr('fill', ({ data }) =>
      data.name && attributes[data.name]
        ? fill.domain([csMin, csMax])(attributes[data.name].Cosine_similarity)
        : stroke
    )
    .attr('r', ({ data }) =>
      data.name && attributes[data.name]
        ? r.domain([mutationMin, mutationMax])(attributes[data.name].Mutations)
        : null
    )
    .on('mouseover', mouseover)
    .on('mousemove', mousemove)
    .on('mouseleave', mouseleave);
  // .call(drag(simulation));

  // simulation.on('tick', () => {
  //   link
  //     .attr('x1', (d) => d.source.x)
  //     .attr('y1', (d) => d.source.y)
  //     .attr('x2', (d) => d.target.x)
  //     .attr('y2', (d) => d.target.y);
  //   node.attr('cx', (d) => d.x).attr('cy', (d) => d.y);
  // });

  const csLegend = Legend(fill.domain([csMin, csMax]), {
    title: 'Cosine Similarity',
    // marginLeft: '1rem',
    // marginTop: '1rem',
  });

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
    // console.log('event', e);
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

  return [container.node(), csLegend];
}

// Copyright 2021, Observable Inc.
// Released under the ISC license.
// https://observablehq.com/@d3/color-legend
function Legend(
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
