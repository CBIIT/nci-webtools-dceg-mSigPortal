import { useRef, useEffect } from 'react';
import * as d3 from 'd3';
import {
  plotStyle,
  getTreeFit,
  getOverlayLayout,
  createColorScale,
  createNodeStyles,
  getTooltipHtml,
} from './treeLeaf.utils';

/**
 * @param {object} props
 * @param {{nodes: object[], links: number[][]} | null} props.layout precomputed positions; links are [sourceIndex, targetIndex]
 * @param {object} props.attributes leaf attributes keyed by sample
 */
export default function D3TreeLeaf({
  id = 'treeleaf-plot',
  width = 1000,
  height = 1000,
  onSelect,
  layout,
  attributes,
  form,
  plotTitle,
}) {
  const plotRef = useRef(null);
  const plotHandleRef = useRef(null);
  const { color, searchSamples } = form;

  useEffect(() => {
    if (!plotRef.current) return;
    if (!layout || !attributes) {
      plotRef.current.replaceChildren();
      plotHandleRef.current = null;
      return;
    }
    const plot = createForceDirectedTree(
      { attributes, nodes: layout.nodes, links: layout.links },
      { id, width, height, radius: Math.min(width, height) / 2 },
      { onClick: onSelect }
    );
    plotRef.current.replaceChildren(plot.node);
    plotHandleRef.current = plot;
  }, [layout, attributes, id, width, height, onSelect]);

  // restyles in place so changing color or search keeps the DOM and zoom
  useEffect(() => {
    plotHandleRef.current?.updateStyle({
      form: { color, searchSamples },
      plotTitle,
    });
  }, [
    layout,
    attributes,
    id,
    width,
    height,
    onSelect,
    color,
    searchSamples,
    plotTitle,
  ]);

  return <div ref={plotRef} />;
}

// Copyright 2022 Observable, Inc.
// Released under the ISC license.
// https://observablehq.com/@d3/force-directed-tree
function createForceDirectedTree(
  { attributes, nodes, links },
  {
    id,
    width = 640, // outer width, in pixels
    height = 400, // outer height, in pixels
    margin = 0, // shorthand for margins
    marginTop = margin, // top margin, in pixels
    marginLeft = margin, // left margin, in pixels
    radius, // outer radius
    strokeOpacity = 1, // stroke opacity for links
    strokeLinejoin, // stroke line join for links
    strokeLinecap, // stroke line cap for links
  },
  { onClick }
) {
  const { viewBoxScale, stroke, strokeWidth } = plotStyle;
  const { cx, cy, treeScale } = getTreeFit(nodes, width, height);

  const zoom = d3.zoom().on('zoom', zoomed);
  const container = d3.create('div').style('position', 'relative');
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

  // add zoom reset control
  const zoomReset = container
    .append('button')
    .attr('id', 'treeleaf-zoom-reset')
    .attr('class', 'btn btn-outline-secondary btn-sm')
    .attr('style', 'position: absolute; top: 10px; left: 10px;')
    .text('Reset Zoom')
    .on('click', () =>
      svg.transition().duration(250).call(zoom.transform, d3.zoomIdentity)
    );

  // add tree container
  const treeZoomContainer = svg
    .append('g')
    .attr('id', 'treeleaf-zoom-container');

  // Center the bounding box at the viewBox origin, then scale to fill.
  const treeGroup = treeZoomContainer
    .append('g')
    .attr('id', 'treeleaf-tree')
    .attr(
      'transform',
      `translate(${-cx * treeScale}, ${-cy * treeScale}) scale(${treeScale})`
    );

  // add lines
  const link = treeGroup
    .append('g')
    .attr('id', 'treeleaf-links')
    .attr('fill', 'none')
    .attr('stroke', stroke)
    .attr('stroke-opacity', strokeOpacity)
    .attr('stroke-linecap', strokeLinecap)
    .attr('stroke-linejoin', strokeLinejoin)
    .attr('stroke-width', strokeWidth * 2)
    .selectAll('path')
    .data(links)
    .join('line')
    .attr('x1', ([source]) => nodes[source].x)
    .attr('y1', ([source]) => nodes[source].y)
    .attr('x2', ([, target]) => nodes[target].x)
    .attr('y2', ([, target]) => nodes[target].y);

  // one listener per event on the group instead of per circle
  const nodeGroup = treeGroup
    .append('g')
    .attr('id', 'treeleaf-nodes')
    .attr('fill', '#fff')
    .attr('stroke', '#000')
    .attr('stroke-width', 1.5)
    .on('mouseover', mouseover)
    .on('mousemove', mousemove)
    .on('mouseout', mouseleave)
    .on('click', click);

  // r=0 circles are never painted, so only visible leaves get DOM elements
  const circles = nodeGroup
    .selectAll('circle')
    .data(nodes.filter((d) => d.r > 0))
    .join('circle')
    .attr('id', (d) => d.name)
    .attr('stroke', stroke)
    .attr('stroke-width', strokeWidth)
    .attr('r', (d) => d.r)
    .attr('cx', (d) => d.x)
    .attr('cy', (d) => d.y)
    .attr('class', 'c-pointer');

  function zoomed({ transform }) {
    treeZoomContainer.attr('transform', transform);
  }

  // add tooltips
  const tooltip = container
    .append('div')
    .style('position', 'absolute')
    .style('visibility', 'hidden')
    .style('pointer-events', 'none')
    .attr('class', 'bg-light border rounded p-1');

  function mouseover() {
    tooltip.style('visibility', 'visible');
  }

  function mouseleave() {
    tooltip.style('visibility', 'hidden');
  }

  function mousemove(e) {
    const sample = d3.select(e.target).datum().name;
    const [px, py] = d3.pointer(e, container.node());
    tooltip
      .html(getTooltipHtml(sample, attributes[sample]))
      .style('left', px + 12 + 'px')
      .style('top', py + 12 + 'px');
  }

  function click(e) {
    const sample = d3.select(e.target).datum().name;
    const data = attributes[sample];
    if (typeof onClick === 'function') {
      onClick(data);
    }
  }

  // append title
  const overlay = getOverlayLayout(width, height);
  const titleText = svg
    .append('g')
    .attr('id', 'treeleaf-title')
    .attr('transform', `translate(0, ${overlay.titleY})`)
    .append('text')
    .attr('x', 0)
    .attr('y', 20)
    .attr('fill', 'black')
    .attr('text-anchor', 'middle')
    .style('font-weight', 'bold')
    .style('font-size', '24px');

  function updateStyle({ form, plotTitle }) {
    const colorFill = createColorScale(attributes, form.color);
    const { getColor, getOpacity } = createNodeStyles({
      attributes,
      form,
      colorScale: colorFill,
    });

    circles
      .attr('fill', (d) => getColor(d.name))
      .attr('opacity', (d) => getOpacity(d.name));

    titleText.text(plotTitle);

    svg.selectAll('#treeleaf-legend, defs').remove();
    const legendParams = {
      svg,
      color: colorFill,
      title: form.color.label,
      x: overlay.legendX,
      y: overlay.legendY,
    };
    form.color.continuous
      ? continuousLegend(legendParams)
      : categoricalLegend(legendParams);
  }

  return { node: container.node(), updateStyle };
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
  size = 25,
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
    .style('font-size', '15px')
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
