import React from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import Description from '../../controls/description/description';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { actions as exposureActions } from '../../../services/store/exposure';
import { actions as modalActions } from '../../../services/store/modal';
import Plot from '../../controls/plot/plot';
import CustomSelect from '../../controls/select/select';
import Debug from '../../controls/debug/debug';
import { NavHashLink } from 'react-router-hash-link';

const actions = { ...exposureActions, ...modalActions };
const { Group } = Form;

export default function MSIndividual({ calculateIndividual }) {
  const dispatch = useDispatch();
  const exposure = useSelector((state) => state.exposure);
  const { sample, plotPath, debugR, err, loading } = exposure.msIndividual;
  const { projectID, publicSampleOptions, userSampleOptions, source } =
    exposure.exposureState;

  const mergeMsIndividual = (state) =>
    dispatch(actions.mergeExposure({ msIndividual: state }));

  return (
    <div>
      <div className="p-3">
        <b>Mutational Signature in Individual Sample</b>
        <Description
          less="The following plot is used to visualize the signature decomposition in individual samples. Select the [Sample Name] and click the [Recalculate] button to
          visualize the signature deconvolution of the selected sample."
          more={
            <>
              <p className="mt-3">
                The combination plot shows the original mutational profile, the
                deconstructed mutational profile, and the difference of each
                mutation type between these two profiles, mutational signature
                profiles and proportion of each contributed signature detected
                in the selected sample. Two measurements (RSS and Cosine
                Similarity) of evaluating signature deconvolution are shown on
                the top of this plot.¸
              </p>
              <p>
                Residual Sum of Squares (RSS) measures the discrepancy between
                two mutational profiles. Cosine similarity measures how similar
                two mutational profiles are. For example, two identical
                mutational profiles will have RSS = 0 and Cosine similarity = 1.
                For additional information about RSS and cosine similarity,
                click{' '}
                <NavHashLink to="/faq#cosine-similarity">here</NavHashLink>.
              </p>
            </>
          }
        />
      </div>
      <hr />
      <Form className="p-3">
        <LoadingOverlay active={loading} />
        <Row className="">
          <Col lg="auto">
            <Group>
              <CustomSelect
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
          <Col lg="auto" className="d-flex">
            <Button
              className="mt-auto mb-3"
              disabled={!sample || (source == 'user' && !projectID)}
              variant="primary"
              onClick={calculateIndividual}
            >
              Recalculate
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
              title="Mutational Signature in Individual Samples"
              downloadName={plotPath.split('/').slice(-1)[0]}
              plotPath={`api/results/${plotPath}`}
              height="900px"
            />
            <p className="p-3">
              The combination plot shows the original mutational profile, the
              deconstructed mutational profile, and the difference of each
              mutation type between these two profiles, mutational signature
              profiles and proportion of each contributed signature detected in
              the selected sample. Two measurements (Residual Sum of Squares,
              RSS, and Cosine Similarity) for evaluating the signature
              deconvolution are shown on the top of this plot.
            </p>
            <p className="p-3">
              RSS measures the discrepancy between two mutational profiles.
              Cosine similarity measures how similar two mutational profiles
              are. For example, two identical mutational profiles will have RSS
              = 0 and Cosine similarity = 1. For additional information about
              RSS and cosine similarity, click{' '}
              <NavHashLink to="/faq#cosine-similarity">here</NavHashLink>.
            </p>
          </>
        )}
        {/* <Debug msg={debugR} /> */}
      </div>
    </div>
  );
}
