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
import { actions as exposureActions } from '../../../../services/store/exploration';
import { useMsIndividualQuery } from './apiSlice';
import Select from '../../../controls/select/selectForm';
import Plotly from '../../../controls/plotly/plot/plot';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { NavHashLink } from 'react-router-hash-link';
import MsIndividualForm from './msIndividual-form';

const { Group, Check } = Form;
const actions = { ...exposureActions };

export default function MsIndividualPlot({ state, form, setForm }) {
  const [params, setParams] = useState('');
  const { data, error, isFetching } = useMsIndividualQuery(params, {
    skip: !params,
  });

  const { sample } = form;
  const { study, strategy, signatureSetName, cancer, useAllCancer, id } = state;

  useEffect(() => {
    let params_activity, params_signature, params_spectrum;
    if (sample && id) {
      params_activity = {
        sample: sample.value,
        userId: id,
      };
      params_signature = {
        userId: id,
      };
      params_spectrum = {
        sample: sample.value,
        userId: id,
      };
      setParams({
        params_activity,
        params_signature,
        params_spectrum,
      });
    } else if (sample && study) {
      params_activity = {
        sample: sample.value,
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        ...(!useAllCancer && { cancer: cancer.value }),
      };
      params_signature = {
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
      };
      params_spectrum = {
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
      <LoadingOverlay active={isFetching} />
      {error && <p className="m-3 alert alert-warning">{error}</p>}
      {data && data.error && (
        <p className="m-3 alert alert-warning">{data.error}</p>
      )}
      {data && !data.error && (
        <div>
          <Plotly
            className="w-100"
            data={data.traces}
            layout={data.layout}
            config={data.config}
          />
          <br />
          <div>
            {' '}
            <p>
              The combination plot shows the original mutational profile, the
              deconstructed mutational profile, and the difference of each
              mutation type between these two profiles, mutational signature
              profiles and proportion of each contributed signature detected in
              the selected sample. Two measurements (Residual Sum of Squares,
              RSS, and Cosine Similarity) for evaluating the signature
              deconvolution are shown on the top of this plot.
            </p>
            <p>
              RSS measures the discrepancy between two mutational profiles.
              Cosine similarity measures how similar two mutational profiles
              are. For example, two identical mutational profiles will have RSS
              = 0 and Cosine similarity = 1. For additional information about
              RSS and cosine similarity, click{' '}
              <NavHashLink to="/faq#cosine-similarity">here</NavHashLink>.
            </p>
          </div>
        </div>
      )}
    </div>
  );
}
