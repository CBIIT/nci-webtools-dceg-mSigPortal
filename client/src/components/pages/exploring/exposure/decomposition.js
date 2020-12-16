import React, { useEffect } from 'react';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchExpDecomposition,
} from '../../../../services/store';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';

export default function Decomposition() {
  const { plotPath, plotURL, txtPath, debugR, err, loading } = useSelector(
    (state) => state.expDecomposition
  );
  const { projectID } = useSelector((state) => state.visualizeResults);

  useEffect(() => {
    if (plotPath) setRPlot(plotPath);
    else clearPlot();
  }, [plotPath]);

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
    if (plotURL) URL.revokeObjectURL(plotURL);
    dispatchExpDecomposition({ plotPath: '', plotURL: '' });
  }

  return (
    <div>
      <LoadingOverlay active={loading} />
      {!err && !plotURL && <p>Please calculate using the left side panel.</p>}
      {err && (
        <div>
          <p>An error has occured. Check the debug section for more info.</p>
          <p>Error: {err}</p>
        </div>
      )}
      {plotURL && (
        <Plot
          plotName={plotPath.split('/').slice(-1)[0]}
          plotURL={plotURL}
          txtPath={projectID + txtPath}
        />
      )}
      <Debug msg={debugR} />
    </div>
  );
}
