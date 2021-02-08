import React, { useEffect } from 'react';
import { useSelector, useDispatch } from 'react-redux';
import { actions as exploringActions } from '../../../../services/store/exploring';
import { actions as modalActions } from '../../../../services/store/modal';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';

const actions = { ...exploringActions, ...modalActions };

export default function MsDecomposition() {
  const dispatch = useDispatch();
  const exploring = useSelector((state) => state.exploring);
  const {
    plotPath,
    plotURL,
    txtPath,
    debugR,
    err,
    loading,
  } = exploring.msDecomposition;
  const { projectID } = exploring.exposure;
  const mergeExploring = (state) =>
    dispatch(actions.mergeExploring({ exploring: state }));
  const mergeMsDecomposition = (state) =>
    dispatch(actions.mergeExploring({ msDecomposition: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

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
          mergeMsDecomposition({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        mergeError(err.message);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      mergeMsDecomposition({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      mergeMsDecomposition({ plotURL: '' });
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
            title="Evaluating the Performance of Mutational Signature Decomposition"
            downloadName={plotPath.split('/').slice(-1)[0]}
            plotURL={plotURL}
            txtPath={projectID + txtPath}
          />
        </>
      )}
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
