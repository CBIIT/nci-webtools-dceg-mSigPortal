import React, { useEffect } from 'react';
import { useSelector, useDispatch } from 'react-redux';
import { actions as exploringActions } from '../../../../services/store/exploring';
import { actions as modalActions } from '../../../../services/store/modal';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';

const actions = { ...exploringActions, ...modalActions };

export default function TmbSignatures() {
  const dispatch = useDispatch();
  const exploring = useSelector((state) => state.exploring);
  const { plotPath, plotURL, debugR, err } = exploring.tmbSignatures;
  const { projectID } = exploring.exposure;
  const mergeExploring = (state) =>
    dispatch(actions.mergeExploring({ exploring: state }));
  const mergeTmbSignatures = (state) =>
    dispatch(actions.mergeExploring({ tmbSignatures: state }));
  const mergeError = (state) => dispatch(actions.mergeModal({ error: state }));

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
          mergeTmbSignatures({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        mergeError({ visible: true, message: err.message });
;
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      mergeTmbSignatures({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      mergeTmbSignatures({ plotURL: '' });
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
