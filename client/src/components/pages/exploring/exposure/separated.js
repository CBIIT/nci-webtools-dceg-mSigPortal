import React, { useEffect } from 'react';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchExpSeparated,
} from '../../../../services/store';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';

export default function Separated() {
  const { plotPath, plotURL, debugR, err } = useSelector(
    (state) => state.expSeparated
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

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      dispatchExpSeparated({ plotPath: '', plotURL: '' });
    }
  }

  return (
    <div>
      <div className="p-3">
        <p>Tumor Mutational Burden Separated by Signatures</p>
        {!err && !plotPath && (
          <p>Please calculate using the left side panel.</p>
        )}
        {err && <p className="text-danger">{err}</p>}
      </div>
      {plotPath && (
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
