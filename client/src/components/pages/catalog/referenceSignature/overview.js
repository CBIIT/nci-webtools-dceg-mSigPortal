import React, { useEffect, useState } from 'react';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Description from '../../../controls/description/description';
import SvgContainer from '../../../controls/svgContainer/svgContainer';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../../services/store/catalog';
import { actions as modalActions } from '../../../../services/store/modal';

const actions = { ...catalogActions, ...modalActions };

export default function Overview() {
  const dispatch = useDispatch();
  const mergeState = async (state) =>
    dispatch(actions.mergeCatalog({ sigRefSig: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));
  const store = useSelector((state) => state.catalog);
  const { plotPath } = store.sigRefSig;

  const [loading, setLoading] = useState(false);

  useEffect(() => {
    async function getPlot() {
      setLoading(true);
      try {
        const plotBlob = await (
          await fetch(`web/getImageS3`, {
            method: 'POST',
            headers: {
              Accept: 'image/svg',
              'Content-Type': 'application/json',
            },
            body: JSON.stringify({
              path: `msigportal/Database/Others/mSigPortalReferenceSignatures.svg`,
            }),
          })
        ).blob();

        const plotURL = URL.createObjectURL(plotBlob);
        mergeState({ plotPath: plotURL });
      } catch (err) {
        mergeState({ plotPath: err });
      }
      setLoading(false);
    }
    if (!plotPath) getPlot();
  }, [plotPath]);

  return (
    <div id="rsPlot">
      <LoadingOverlay active={loading} />
      <div className="p-3">
        The pie charts below display the current reference signatures (RS)
        available in mSigPortal for both human (GRCh37/38) and mouse genome
        (GRCm38). Each pie chart represents a given mutational signature defined
        by the profile type (SBS, DBS, ID, RS) and its respective matrix size.
        Each signature set included in mSigPortal is denoted by a color in the
        legend on the right. The numbers and coloring in each chart represent
        the number of signatures included and the signature source,
        respectively.
      </div>
      <hr />
      <SvgContainer
        className="p-3"
        downloadName={plotPath.split('/').slice(-1)[0]}
        plotPath={plotPath}
        height="800px"
      />
    </div>
  );
}
