import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { actions as exposureActions } from '../../../services/store/exposure';
import { actions as modalActions } from '../../../services/store/modal';
import Plot from '../../controls/plot/plot';
import Select from '../../controls/select/select';
import Debug from '../../controls/debug/debug';

const actions = { ...exposureActions, ...modalActions };
const { Group } = Form;

export default function MSIndividual({ calculateIndividual }) {
  const dispatch = useDispatch();
  const exposure = useSelector((state) => state.exposure);
  const { sample, plotPath, debugR, err, loading } = exposure.msIndividual;
  const {
    projectID,
    publicSampleOptions,
    userSampleOptions,
    source,
  } = exposure.exposureState;
  const mergeExposure = (state) =>
    dispatch(actions.mergeExposure({ exposure: state }));
  const mergeMsIndividual = (state) =>
    dispatch(actions.mergeExposure({ msIndividual: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  // choose inital sample
  useEffect(() => {
    if (source == 'public') {
      mergeMsIndividual({ sample: publicSampleOptions[0] });
    } else {
      mergeMsIndividual({ sample: userSampleOptions[0] || '' });
    }
  }, [publicSampleOptions, userSampleOptions, source]);

  return (
    <div>
      <div className="p-3">
        <b>Mutational Signature in Individual Sample</b>
        <p className="m-0">
          This page allows you to visualize the signature decomposition in
          individual samples. Selected the sample name and click “Calculate”
          button to visualize the new sample. In this plot, it will show the
          original mutational profile, deconstructed mutational profile and the
          difference between these two profiles. Also at the top of the plot are
          measurements for RSS and cosine similarity. RSS is the Residual Sum of
          Squares. It measures the discrepancy between two profiles. Cosine
          similarity is how similar the mutational profiles are to one another.
          For additional information about RSS and cosine similarity, click
          here. In addition, all relevant signature profiles from signature
          decomposition will present on the bottom. A simple formula (on the
          bottom) as well as a bar plot (on the left) will also show the
          signature contribution.
        </p>
      </div>
      <hr />
      <Form className="p-3">
        <LoadingOverlay active={loading} />
        <Row className="">
          <Col lg="3">
            <Group>
              <Select
                disabled={
                  source == 'public'
                    ? !publicSampleOptions.length
                    : !userSampleOptions.length
                }
                id="msIndSample"
                label="Sample Name"
                value={sample}
                options={
                  source == 'public' ? publicSampleOptions : userSampleOptions
                }
                onChange={(name) => mergeMsIndividual({ sample: name })}
              />
            </Group>
          </Col>
          <Col />
          <Col lg="2" className="d-flex">
            <Button
              className="ml-auto mb-auto"
              disabled={!sample || (source == 'user' && !projectID)}
              variant="primary"
              onClick={calculateIndividual}
            >
              Calculate
            </Button>
          </Col>
        </Row>
      </Form>
      <div id="msIndividualPlot">
        {err && (
          <div>
            <hr />
            <p className="p-3 text-danger">{err}</p>
          </div>
        )}
        {plotPath && (
          <>
            <hr />
            <Plot
              className="p-3"
              title="Mutational Signature in Individual Sample"
              downloadName={plotPath.split('/').slice(-1)[0]}
              plotPath={`api/results/${projectID}${plotPath}`}
              maxHeight="700px"
            />
          </>
        )}
        {/* <Debug msg={debugR} /> */}
      </div>
    </div>
  );
}
