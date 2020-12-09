import React, { useEffect } from 'react';
import { useSelector } from 'react-redux';
import { dispatchError, dispatchExpTumor } from '../../../../services/store';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';

export default function Tumor() {
  const { plotPath, plotURL, debugR, err } = useSelector(
    (state) => state.expTumor
  );

  useEffect(() => {
    if (plotPath) setRPlot(plotPath);
    else clearPlot();
  }, [plotPath]);

  async function setRPlot(plotPath) {
    if (plotPath) {
      try {
        const response = await fetch(`api/results/${plotPath}`);
        if (!response.ok) {
          // console.log(await response.json());
        } else {
          const pic = await response.blob();
          const objectURL = URL.createObjectURL(pic);

          if (plotURL) URL.revokeObjectURL(plotURL);
          dispatchExpTumor({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchExpTumor({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) URL.revokeObjectURL(plotURL);
    dispatchExpTumor({ plotPath: '', plotURL: '' });
  }

  return (
    <div>
      {!err && !plotURL && <p>Please calculate using the left side panel.</p>}
      {err && <p>{err}</p>}
      {plotURL && (
        <Plot plotName={plotPath.split('/').slice(-1)[0]} plotURL={plotURL} />
      )}
      <Debug msg={debugR} />
    </div>
  );
}
