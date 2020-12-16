import React, { useEffect } from 'react';
import { Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { dispatchError, dispatchExpRefSig } from '../../../../services/store';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';

export default function ReferenceSignatures({ submitR }) {
  const { plotPath, plotURL, debugR, err, loading } = useSelector(
    (state) => state.expRefSig
  );
  // const { displayTab } = useSelector((state) => state.exploring);

  // useEffect(() => {
  //   if (!loading && !plotPath && displayTab == 'signatureExploring') {
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
  //     dispatchError(err);
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
  //       dispatchError(err);
  //     }
  //   } else {
  //     if (plotURL) URL.revokeObjectURL(plotURL);
  //     dispatchExpRefSig({ err: true, plotURL: '' });
  //   }
  // }

  return (
    <div>
      <div id="plot">
        <div style={{ display: err ? 'block' : 'none' }}>
          <p>An error has occured. Check the debug section for more info.</p>
        </div>
        <div style={{ display: plotURL ? 'block' : 'none' }}>
          <Plot
            plotName={plotURL.split('/').slice(-1)[0]}
            plotURL={plotURL}
            maxHeight="600px"
          />
        </div>
      </div>
    </div>
  );
}
