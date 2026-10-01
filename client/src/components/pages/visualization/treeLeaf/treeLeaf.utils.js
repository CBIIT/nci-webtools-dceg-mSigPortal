import * as d3 from 'd3';

function downloadBlob(blob, filename) {
  const url = URL.createObjectURL(blob);
  const downloadLink = document.createElement('a');
  downloadLink.href = url;
  downloadLink.download = filename;
  document.body.appendChild(downloadLink);
  downloadLink.click();
  document.body.removeChild(downloadLink);
  setTimeout(() => URL.revokeObjectURL(url), 0);
}

export function exportSvg(selector, filename) {
  const svgString = new XMLSerializer().serializeToString(
    document.querySelector(selector)
  );
  downloadBlob(new Blob([svgString], { type: 'image/svg+xml' }), filename);
}

export function groupBy(array, key) {
  const result = Object.create(null);
  for (const item of array) {
    result[item[key]] = item;
  }
  return result;
}

export const plotStyle = {
  viewBoxScale: 1.3,
  fillFactor: 0.7,
  stroke: '#666',
  strokeWidth: 0.3,
  highlightColor: 'yellow',
};

export function getPlotTitle({ isUser, form, studyLabel, signatureSetName }) {
  let plotTitle = isUser
    ? `User Data - ${form.color.label}`
    : `${studyLabel} - ${form.color.label}`;

  if (
    !isUser &&
    form.color.label === 'Dominant Signature' &&
    signatureSetName
  ) {
    plotTitle += ' - ' + signatureSetName;
  }
  return plotTitle;
}

/**
 * Fits the node bounding box (including circle radii, so few-node plots don't
 * clip) into the viewBox with a margin.
 */
export function getTreeFit(nodes, width, height) {
  const { viewBoxScale, fillFactor } = plotStyle;
  const xMin = d3.min(nodes, (d) => d.x - (d.r || 0));
  const xMax = d3.max(nodes, (d) => d.x + (d.r || 0));
  const yMin = d3.min(nodes, (d) => d.y - (d.r || 0));
  const yMax = d3.max(nodes, (d) => d.y + (d.r || 0));
  const bboxW = Math.max(xMax - xMin, 1);
  const bboxH = Math.max(yMax - yMin, 1);
  const cx = (xMin + xMax) / 2;
  const cy = (yMin + yMax) / 2;
  const treeScale =
    fillFactor *
    Math.min((width * viewBoxScale) / bboxW, (height * viewBoxScale) / bboxH);
  return { cx, cy, treeScale };
}

/** Title and legend anchors in viewBox coordinates. */
export function getOverlayLayout(width, height) {
  const { viewBoxScale } = plotStyle;
  return {
    titleY: (-height / 2 + 10) * viewBoxScale,
    legendX: (width / 2 - 40) * viewBoxScale,
    legendY: (-height / 2 + 20) * viewBoxScale,
  };
}

export function createColorScale(attributes, color) {
  const colorValues = Object.values(attributes).map((e) => e[color.value]);
  return color.continuous
    ? d3
        .scaleSequential(d3.interpolateRgb('white', 'steelblue'))
        .domain([d3.min(colorValues), d3.max(colorValues)])
    : d3.scaleOrdinal(d3.schemeCategory10).domain(colorValues);
}

export function createNodeStyles({ attributes, form, colorScale }) {
  const { stroke, highlightColor } = plotStyle;
  const searchValues = form.searchSamples?.map((s) => s.value) || [];

  function getColor(name) {
    if (name && searchValues.includes(name)) {
      return highlightColor;
    }

    return name && attributes[name]
      ? colorScale(attributes[name][form.color.value])
      : stroke;
  }

  function getOpacity(name) {
    if (!searchValues?.length || (name && searchValues.includes(name))) {
      return 1;
    }

    return 0.8;
  }

  return { getColor, getOpacity };
}

function escapeHtml(value) {
  return String(value)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&#39;');
}

export function getTooltipHtml(sample, data = {}) {
  const cosine =
    data.Cosine_similarity != null && !Number.isNaN(+data.Cosine_similarity)
      ? (+data.Cosine_similarity).toFixed(3)
      : 'Unavailable';
  const rows = [
    ['Sample', sample],
    ['Cancer Type', data.Cancer_Type],
    ['Cosine Similarity', cosine],
    ['Dominant Mutation', data.Dmut],
    ['Mutations', data.Mutations],
  ];
  return `<div class="text-start">${rows
    .map(
      ([label, value]) =>
        `<div>${label}: ${escapeHtml(value ?? 'Unavailable')}</div>`
    )
    .join('')}</div>`;
}
