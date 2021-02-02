import React, { useEffect } from 'react';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchTmbSignatures,
} from '../../../../services/store';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';

export default function TmbSignatures() {
  const { plotPath, plotURL, debugR, err } = useSelector(
    (state) => state.tmbSignatures
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
          dispatchTmbSignatures({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchTmbSignatures({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      dispatchTmbSignatures({ plotURL: '' });
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
            plotName="Tumor Mutational Burden Separated by Signatures"
            plotURL={plotURL}
          />
        </>
      )}
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
