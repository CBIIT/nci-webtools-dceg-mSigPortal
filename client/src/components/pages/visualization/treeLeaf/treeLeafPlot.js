import { useRef, useEffect } from 'react';
import { forceSimulation, forceLink } from 'd3-force';
import * as d3 from 'd3';
import cloneDeep from 'lodash/cloneDeep';
import { useRecoilValue } from 'recoil';
import { treeLeafData } from './treeLeaf.state';

export default function D3TreeLeaf({ width = 800, height = 500, ...props }) {
  const { links, nodes } = useRecoilValue(treeLeafData);
  const nodeRef = useRef(null);

  useEffect(() => {
    if (nodeRef.current && nodes.length && links.length) {
      const treeLeaf = createTreeLeaf(
        cloneDeep(nodes),
        cloneDeep(links),
        width,
        height
      );
      nodeRef.current.replaceChildren(treeLeaf);
    }
  }, [links, nodes, nodeRef]);

  return <div ref={nodeRef} {...props} />;
}

export function createTreeLeaf({ nodes, links, width, height }) {
  const simulation = forceSimulation(nodes)
    .force('link', d3.forceLink(links))
    .force('charge', d3.forceManyBody())
    .force('x', d3.forceX())
    .force('y', d3.forceY());
  // .force('center', d3.forceCenter());

  const svg = d3
    .create('svg')
    .attr('viewBox', [-width / 2, -height / 2, width, height])
    .attr('width', width)
    .attr('height', height);

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
    .attr('r', 3.5)
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
