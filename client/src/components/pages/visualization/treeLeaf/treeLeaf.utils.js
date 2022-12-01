export function exportSvg(selector, filename) {
  const svgString = new XMLSerializer().serializeToString(document.querySelector(selector));
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
