import React, { useEffect } from 'react';
import { useSelector } from 'react-redux';
import { dispatchError, dispatchExpTumor } from '../../../../services/store';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';

export default function Tumor() {
  const { plotPath, plotURL, debugR, err } = useSelector(
    (state) => state.expTumor
  );
  const { projectID } = useSelector((state) => state.expExposure);

  useEffect(() => {
    if (plotPath) setRPlot(plotPath);
    else clearPlot();
  }, [plotPath, err, debugR, projectID]);

  async function setRPlot(plotPath) {
    if (plotPath) {
      try {
        const response = await fetch(`api/results/${projectID}${plotPath}`);
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
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      dispatchExpTumor({ plotPath: '', plotURL: '' });
    }
  }

  return (
    <div>
      <div className="p-3">
        <p>Tumor Mutational Burden</p>
        {!err && !plotPath && (
          <p>Please calculate using the left side panel.</p>
        )}
        {err && <p className="text-danger">{err}</p>}
      </div>
      {plotURL && (
        <>
          <hr />
          <Plot
            className="p-3"
            plotName={plotPath.split('/').slice(-1)[0]}
            plotURL={plotURL}
          />
        </>
      )}
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
