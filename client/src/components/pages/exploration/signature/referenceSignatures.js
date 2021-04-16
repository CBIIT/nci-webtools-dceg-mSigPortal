import React, { useEffect } from 'react';
import { Row, Col, Button } from 'react-bootstrap';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import { useSelector, useDispatch } from 'react-redux';
import { actions as explorationActions } from '../../../../services/store/exploration';
import { actions as modalActions } from '../../../../services/store/modal';

const actions = { ...explorationActions, ...modalActions };

export default function ReferenceSignatures({ submitR }) {
  const dispatch = useDispatch();
  const exploration = useSelector((state) => state.exploration);
  const { plotPath, plotURL, debugR, err, loading } = exploration.sigRefSig;
  // const { displayTab } = useSelector((state) => state.exploration);

  // useEffect(() => {
  //   if (!loading && !plotPath && displayTab == 'signatureExploration') {
  //     calculateR('referenceSignatures', {});
  //   }
  // }, [plotPath, displayTab]);

  // async function calculateR(fn, args) {
  //   dispatchExpRefSig({
  //     loading: true,
  //     err: false,
  //     debugR: '',
  //   });

  //   try {
  //     const response = await submitR(fn, args);
  //     if (!response.ok) {
  //       const err = await response.json();

  //       dispatchExpRefSig({
  //         loading: false,
  //         debugR: err,
  //       });
  //     } else {
  //       const { debugR, output } = await response.json();

  //       dispatchExpRefSig({
  //         debugR: debugR,
  //         loading: false,
  //         plotPath: output.plotPath,
  //       });
  //       setRPlot(output.plotPath, 'within');
  //     }
  //   } catch (err) {
  //     mergeError(err.message );
  //     dispatchExpRefSig({ loading: false });
  //   }
  // }

  // async function setRPlot(plotPath) {
  //   if (plotPath) {
  //     try {
  //       const response = await fetch(`api/results/${projectID}${plotPath}`);
  //       if (!response.ok) {
  //         // console.log(await response.json());
  //       } else {
  //         const pic = await response.blob();
  //         const objectURL = URL.createObjectURL(pic);

  //         if (plotURL) URL.revokeObjectURL(plotURL);
  //         dispatchExpRefSig({
  //           plotURL: objectURL,
  //         });
  //       }
  //     } catch (err) {
  //       mergeError(err.message );
  //     }
  //   } else {
  //     if (plotURL) URL.revokeObjectURL(plotURL);
  //     dispatchExpRefSig({ err: true, plotURL: '' });
  //   }
  // }

  return (
    <div id="plot">
      <div style={{ display: err ? 'block' : 'none' }}>
        <p>An error has occured. Please verify your input.</p>
      </div>
      <div style={{ display: plotURL ? 'block' : 'none' }}>
        <p className="p-3">
          The pie charts below displays the reference signatures currently
          available in mSigPortal. Each signature set included in mSigPortal is
          denoted by a color in the legend on the right. Each pie chart
          represents a given mutational signatures defined by specific profile
          type (SBS, DBS, ID, RS) and its respective matrix size. The numbers
          and coloring in each chart represent the number of signatures included
          and the name of the signature source.{' '}
        </p>
        <hr />
        <Plot
          className="p-3"
          downloadName={plotURL.split('/').slice(-1)[0]}
          plotURL={plotURL}
          maxHeight="600px"
        />
      </div>
    </div>
  );
}
