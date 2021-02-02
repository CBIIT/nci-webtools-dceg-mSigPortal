import React, { useEffect } from 'react';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchMsDecomposition,
} from '../../../../services/store';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';

export default function MsDecomposition() {
  const { plotPath, plotURL, txtPath, debugR, err, loading } = useSelector(
    (state) => state.msDecomposition
  );
  const { projectID } = useSelector((state) => state.expExposure);

  useEffect(() => {
    plotPath ? setRPlot(plotPath) : clearPlot();
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
          dispatchMsDecomposition({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchMsDecomposition({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      dispatchMsDecomposition({ plotURL: '' });
    }
  }

  return (
    <div>
      <div>
        {!err && !plotPath && (
          <p className="p-3">Please calculate using the left side panel.</p>
        )}
        {err && (
          <div>
            <hr />
            <p className="p-3 text-danger">{err}</p>
          </div>
        )}
      </div>
      {plotPath && (
        <>
          <Plot
            className="p-3"
            plotName="Evaluating the Performance of Mutational Signature Decomposition"
            plotURL={plotURL}
            txtPath={projectID + txtPath}
          />
        </>
      )}
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
