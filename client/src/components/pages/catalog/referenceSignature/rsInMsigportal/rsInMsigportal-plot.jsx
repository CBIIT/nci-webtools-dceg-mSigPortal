import React, { useState } from 'react';
import { LoadingOverlay } from '../../../../controls/loading-overlay/loading-overlay';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../../../services/store/catalog';
import { actions as modalActions } from '../../../../../services/store/modal';
import Plotly from '../../../../controls/plotly/plot/plot';
import { useRsInMsigportalDataQuery } from './apiSlice';
import { cloneDeep } from 'lodash';

const actions = { ...catalogActions, ...modalActions };

export default function RsInMsigportal() {
  const dispatch = useDispatch();
  const mergeRsInMsigportal = async (state) =>
    dispatch(actions.mergeCatalog({ sigRefSig: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));
  const store = useSelector((state) => state.catalog);
  const { plot } = store.sigRefSig;

  const [params, setParams] = useState(null);

  const [loading, setLoading] = useState(false);

  const { data, error, isFetching } = useRsInMsigportalDataQuery();

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
      {error ? (
        <div className="text-center">
          <div>An error has occured</div>
          <div>{error.message}</div>
        </div>
      ) : (
        data && (
          <Plotly
            className="w-100"
            data={cloneDeep(data.traces)}
            layout={cloneDeep(data.layout)}
            config={cloneDeep(data.config)}
          />
        )
      )}
    </div>
  );
}
