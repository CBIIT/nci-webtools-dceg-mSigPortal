import React from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import Description from '../../controls/description/description';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { actions as exposureActions } from '../../../services/store/exposure';
import { actions as modalActions } from '../../../services/store/modal';
import Plot from '../../controls/plot/plot';
import Select from '../../controls/select/select';

const actions = { ...exposureActions, ...modalActions };

export default function MsBurden({ calculateBurden }) {
  const dispatch = useDispatch();
  const exposure = useSelector((state) => state.exposure);

  const { signatureName, plotPath, debugR, err, loading } = exposure.msBurden;
  const { projectID, signatureNameOptions, userNameOptions, source } =
    exposure.exposureState;

  const mergeMsBurden = (state) =>
    dispatch(actions.mergeExposure({ msBurden: state }));

  return (
    <div>
      <div className="p-3">
        <b>Mutational Signature Burden Across Cancer Types</b>
        <Description
          less="The bar plot below illustrates mutational signature burden across different cancer types with regard to a specific selected signature."
          more="On the y-axis is the number of mutations per Mb (log10) assigned to selected signatures, and the x-axis denotes sample numbers. The green number is the number of samples for each cancer type, and the blue number is the number of samples in that cancer type that detected the selected signature."
        />
      </div>
      <hr />
      <Form className="p-3">
        <LoadingOverlay active={loading} />
        <Row>
          <Col lg="auto">
            <Select
              className="mb-2"
              disabled={source == 'user' && !userNameOptions.length}
              id="acrossSignatureName"
              label="Signature Name"
              value={signatureName}
              options={
                source == 'public' ? signatureNameOptions : userNameOptions
              }
              onChange={(name) => mergeMsBurden({ signatureName: name })}
            />
          </Col>
          <Col lg="auto" className="d-flex">
            <Button
              disabled={!signatureName || (source == 'user' && !projectID)}
              className="mt-auto mb-2"
              variant="primary"
              onClick={calculateBurden}
            >
              Recalculate
            </Button>
          </Col>
        </Row>
      </Form>
      <div id="exposureAcrossPlot">
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
              title="Mutational Signature Burden Across Cancer Types"
              downloadName={plotPath.split('/').slice(-1)[0]}
              plotPath={`api/results/${plotPath}`}
            />
          </>
        )}
        {/* <Debug msg={debugR} /> */}
      </div>
    </div>
  );
}
