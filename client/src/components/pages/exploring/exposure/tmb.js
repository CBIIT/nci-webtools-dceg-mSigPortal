import React, { useEffect } from 'react';
import { useSelector, useDispatch } from 'react-redux';
import { actions as exploringActions } from '../../../../services/store/exploring';
import { actions as modalActions } from '../../../../services/store/modal';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';

const actions = { ...exploringActions, ...modalActions };

export default function TMB() {
  const dispatch = useDispatch();
  const exploring = useSelector((state) => state.exploring);
  const { plotPath, plotURL, debugR, err } = exploring.tmb;
  const { projectID } = exploring.exposure;
  const mergeExploring = (state) =>
    dispatch(actions.mergeExploring({ exploring: state }));
  const mergeTMB = (state) => dispatch(actions.mergeExploring({ tmb: state }));
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
          mergeTMB({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        mergeError({ visible: true, message: err.message });
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      mergeTMB({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      mergeTMB({ plotURL: '' });
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
      {plotURL && (
        <>
          <Plot
            className="p-3"
            plotName="Tumor Mutational Burden"
            plotURL={plotURL}
          />
        </>
      )}
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
