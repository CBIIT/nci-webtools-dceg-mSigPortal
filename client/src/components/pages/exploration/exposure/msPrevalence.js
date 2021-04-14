import React, { useEffect, useState } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { actions as explorationActions } from '../../../../services/store/exploration';
import { actions as modalActions } from '../../../../services/store/modal';
import Plot from '../../../controls/plot/plot';
import Select from '../../../controls/select/select';
import Debug from '../../../controls/debug/debug';

const actions = { ...explorationActions, ...modalActions };
const { Group, Label, Control } = Form;

export default function MsPrevalence({ calculatePrevalence }) {
  const dispatch = useDispatch();
  const exploration = useSelector((state) => state.exploration);
  const {
    mutation,
    plotPath,
    plotURL,
    debugR,
    err,
    loading,
  } = exploration.msPrevalence;
  const { projectID, source } = exploration.exposure;
  const mergeExploration = (state) =>
    dispatch(actions.mergeExploration({ exploration: state }));
  const mergeMsPrevalence = (state) =>
    dispatch(actions.mergeExploration({ msPrevalence: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const [invalidMin, setMin] = useState(false);

  useEffect(() => {
    plotPath ? setRPlot(plotPath) : clearPlot();
  }, [plotPath, err, debugR]);

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
        mergeError(err.message);
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
            <Group
              controlId="prevalenceMutations"
              title="Minimum Number of Mutations Within Each Signature"
            >
              <Label>Minimum Number of Mutations Assigned to Each Signature</Label>
              <Control
                value={mutation}
                placeholder="e.g. 100"
                onChange={(e) => {
                  mergeMsPrevalence({
                    mutation: e.target.value,
                  });
                }}
                isInvalid={invalidMin}
              />
              <Form.Control.Feedback type="invalid">
                Enter a numeric minimum value
              </Form.Control.Feedback>
            </Group>
          </Col>
          <Col />
          <Col lg="2" className="d-flex">
            <Button
              className="ml-auto mb-auto"
              variant="primary"
              onClick={() => {
                if (!mutation || isNaN(mutation)) setMin(true);
                else setMin(false);

                if (mutation && !isNaN(mutation)) calculatePrevalence();
              }}
              disabled={source == 'user' && !projectID}
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
            title="Prevalence of Mutational Signature"
            downloadName={plotPath.split('/').slice(-1)[0]}
            plotURL={plotURL}
          />
        </div>
      </div>
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
