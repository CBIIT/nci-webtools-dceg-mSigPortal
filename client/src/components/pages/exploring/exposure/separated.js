import React, { useEffect } from 'react';
import { useSelector } from 'react-redux';
import { dispatchError, dispatchExpSeparated } from '../../../../services/store';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';

export default function Separated() {
  const { loading } = useSelector((state) => state.expExposure);
  const { plotPath, plotURL, debugR, err } = useSelector(
    (state) => state.expSeparated
  );

  useEffect(() => {
    if (plotPath) setRPlot(plotPath);
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
          dispatchExpSeparated({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchExpSeparated({ err: true, plotURL: '' });
    }
  }

  return (
    <div>
      <LoadingOverlay active={loading} />
      {err && (
        <p>An error has occured. Check the debug section for more info.</p>
      )}
      {plotURL && (
        <Plot plotName={plotPath.split('/').slice(-1)[0]} plotURL={plotURL} />
      )}
      <Debug msg={debugR} />
    </div>
  );
}
