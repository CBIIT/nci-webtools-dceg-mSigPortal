import React, { useEffect } from 'react';
import { Row, Col, Button } from 'react-bootstrap';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Description from '../../../controls/description/description';
import Plot from '../../../controls/plot/plot';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../../services/store/catalog';
import { actions as modalActions } from '../../../../services/store/modal';

const actions = { ...catalogActions, ...modalActions };

export default function ReferenceSignatures({ submitR }) {
  const dispatch = useDispatch();
  const catalog = useSelector((state) => state.catalog);
  const { plotPath, debugR, err, loading } = catalog.sigRefSig;
  // const { displayTab } = useSelector((state) => state.catalog);

  // useEffect(() => {
  //   if (!loading && !plotPath && displayTab == 'signatureCatalog') {
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

  return (
    <div id="plot">
      <div style={{ display: err ? 'block' : 'none' }}>
        <p>An error has occured. Please verify your input.</p>
      </div>
      <div style={{ display: plotPath ? 'block' : 'none' }}>
        <Description
          className="p-3 m-0"
          less="The pie charts below display the current reference signatures (RS) available in mSigPortal for both human (GRCh37/38) and mouse genome (GRCm38)."
          more="Each pie chart represents a given mutational signature defined by profile type (SBS, DBS, ID, RS) and its respective matrix size. Each signature set included in mSigPortal is denoted by a color in the legend on the right. The numbers and coloring in each chart represent the number of signatures included and the signature source, respectively."
        />
        <hr />
        <Plot
          className="p-3"
          downloadName={plotPath.split('/').slice(-1)[0]}
          plotPath={plotPath}
          height="800px"
        />
      </div>
    </div>
  );
}
