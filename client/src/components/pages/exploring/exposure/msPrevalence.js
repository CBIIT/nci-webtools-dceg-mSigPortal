import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { actions as exploringActions } from '../../../../services/store/exploring';
import { actions as modalActions } from '../../../../services/store/modal';
import Plot from '../../../controls/plot/plot';
import Select from '../../../controls/select/select';
import Debug from '../../../controls/debug/debug';

const actions = { ...exploringActions, ...modalActions };
const { Group, Label, Control } = Form;

export default function MsPrevalence({ calculatePrevalence }) {
  const dispatch = useDispatch();
  const exploring = useSelector((state) => state.exploring);
  const {
    mutation, plotPath, plotURL, debugR, err, loading
  } = exploring.msPrevalence;
  const { projectID } = exploring.exposure;
  const mergeExploring = (state) =>
    dispatch(actions.mergeExploring({ exploring: state }));
  const mergeMsPrevalence = (state) =>
    dispatch(actions.mergeExploring({ msPrevalence: state }));
  const mergeError = (state) => dispatch(actions.mergeModal({ error: state }));


  useEffect(() => {
    plotPath ? setRPlot(plotPath) : clearPlot();
  }, [plotPath, err, debugR, projectID]);

  async function setRPlot(plotPath) {
    if (plotPath) {
      try {
        const response = await fetch(`api/results/${projectID}${plotPath}`);
        if (!response.ok) {
          // console.log(await response.json());
        } else {
          const pic = await response.blob();
          const objectURL = URL.createObjectURL(pic);

          if (plotURL) URL.revokeObjectURL(plotURL);
          mergeMsPrevalence({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        mergeError({ visible: true, message: err.message });
;
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      mergeMsPrevalence({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      mergeMsPrevalence({ plotURL: '' });
    }
  }

  return (
    <div>
      <Form noValidate className="p-3">
        <LoadingOverlay active={loading} />
        <Row>
          <Col lg="5">
            <Group controlId="prevalenceMutations">
              <Label>Minimal Number Mutations within in Each Signature</Label>
              <Control
                value={mutation}
                placeholder="e.g. 100"
                onChange={(e) => {
                  mergeMsPrevalence({
                    mutation: e.target.value,
                  });
                }}
                isInvalid={!mutation || isNaN(mutation)}
              />
              <Form.Control.Feedback type="invalid">
                Enter a minimum value
              </Form.Control.Feedback>
            </Group>
          </Col>
          <Col lg="5" />
          <Col lg="2" className="d-flex justify-content-end">
            <Button
              className="mt-auto mb-3"
              variant="primary"
              onClick={calculatePrevalence}
              disabled={!mutation || isNaN(mutation)}
            >
              Calculate
            </Button>
          </Col>
        </Row>
      </Form>
      <div id="exposurePrevalencePlot">
        {err && (
          <div>
            <p className="p-3 text-danger">{err}</p>
          </div>
        )}
        <div style={{ display: plotURL ? 'block' : 'none' }}>
          <hr />
          <Plot
            className="p-3"
            plotName="Prevalence of Mutational Signature"
            plotURL={plotURL}
          />
        </div>
      </div>
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
