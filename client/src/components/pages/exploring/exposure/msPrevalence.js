import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchMsPrevalence,
} from '../../../../services/store';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';

const { Group, Label, Control } = Form;

export default function MsPrevalence({ calculatePrevalence }) {
  const { mutation, plotPath, plotURL, debugR, err, loading } = useSelector(
    (state) => state.msPrevalence
  );
  const { projectID } = useSelector((state) => state.expExposure);

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
          dispatchMsPrevalence({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchMsPrevalence({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      dispatchMsPrevalence({ plotURL: '' });
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
                  dispatchMsPrevalence({
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
