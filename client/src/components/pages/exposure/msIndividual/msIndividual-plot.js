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

import Select from '../../../controls/select/selectForm';
import Plotly from '../../../controls/plotly/plot/plot';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
const { Group, Check } = Form;
const actions = { ...exposureActions };

export default function MsIndividualPlot({ state }) {
  const dispatch = useDispatch();
  const exposure = useSelector((state) => state.exposure);
  const { sample, plotPath, debugR, err, loading } = exposure.msIndividual;
  const { projectID, publicSampleOptions, userSampleOptions, source } =
    exposure.main;

  const mergeMsIndividual = (state) =>
    dispatch(actions.mergeExposure({ msIndividual: state }));
  return <Form className="p-3"></Form>;
}
