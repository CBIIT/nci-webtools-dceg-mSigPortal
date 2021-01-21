import React, { useEffect } from 'react';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchExpDecomposition,
} from '../../../../services/store';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';

export default function Decomposition() {
  const { plotPath, plotURL, txtPath, debugR, err, loading } = useSelector(
    (state) => state.expDecomposition
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
          dispatchExpDecomposition({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchExpDecomposition({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      dispatchExpDecomposition({ plotPath: '', plotURL: '' });
    }
  }

  return (
    <div>
      <div className="p-3">
        <p>Evaluating the Performance of Mutational Signature Decomposition</p>
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
            txtPath={projectID + txtPath}
          />
        </>
      )}
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
