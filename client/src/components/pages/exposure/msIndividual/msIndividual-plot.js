import { useState, useEffect } from 'react';
import {
  Form,
  Row,
  Col,
  Button,
  OverlayTrigger,
  Popover,
} from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { actions as exposureActions } from '../../../../services/store/exposure';
import { useMsIndividualQuery } from './apiSlice';
import Select from '../../../controls/select/selectForm';
import Plotly from '../../../controls/plotly/plot/plot';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
const { Group, Check } = Form;
const actions = { ...exposureActions };

export default function MsIndividualPlot({ state, form }) {
  const [params, setParams] = useState('');
  const { data, error, isFetching } = useMsIndividualQuery(params, {
    skip: !params,
  });

  const { sample } = form;
  const { study, strategy, signatureSetName, cancer, useAllCancer, id } = state;

  useEffect(() => {
    if (sample && id) {
      setParams({
        sample: sample.value,
        userId: id,
      });
    } else if (sample && study) {
      const params_activity = {
        sample: sample.value,
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        ...(!useAllCancer && { cancer: cancer.value }),
      };
      const params_signature = {
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
      };
      const params_spectrum = {
        study: study.value,
        strategy: strategy.value,
        sample: sample.value,
        ...(!useAllCancer && { cancer: cancer.value }),
      };
      setParams({
        params_activity,
        params_signature,
        params_spectrum,
      });
    }
  }, [sample, id]);
  return (
    <div>
      {/* <div className="p-3 text-center">
        <h6>
          <b>Mutational Signature in Individual Samples</b>
        </h6>
      </div> */}
      <LoadingOverlay active={isFetching} />
      {error && <p className="p-3 text-danger">{error}</p>}
      {data && (
        <Plotly
          className="w-100"
          data={data.traces}
          layout={data.layout}
          config={data.config}
        />
      )}
    </div>
  );
}
