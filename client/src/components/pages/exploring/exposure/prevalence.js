import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchExpPrevalence,
} from '../../../../services/store';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';
import Select from '../../../controls/select/select';

const { Group, Label, Control, Text } = Form;

export default function Tumor({ calculatePrevalence }) {
  const { mutation, plotPath, plotURL, debugR, err, loading } = useSelector(
    (state) => state.expPrevalence
  );
  const { projectID } = useSelector((state) => state.visualizeResults);

  useEffect(() => {
    if (plotPath) setRPlot(plotPath);
    else clearPlot();
  }, [plotPath]);

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
          dispatchExpPrevalence({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchExpPrevalence({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) URL.revokeObjectURL(plotURL);
    dispatchExpPrevalence({ plotPath: '', plotURL: '' });
  }

  return (
    <div>
      <Form>
        <LoadingOverlay active={loading} />
        <Row className="">
          <Col sm="4">
            <Group controlId="prevalenceMutations">
              <Label>Minimal Number Mutations within in Each Signature</Label>
              <Control
                value={mutation}
                placeholder="e.g. 100"
                onChange={(e) => {
                  dispatchExpPrevalence({
                    mutation: e.target.value,
                  });
                }}
              />
              {/* <Text className="text-muted">(Ex. NCG>NTG)</Text> */}
            </Group>
          </Col>
          <Col sm="6" />
          <Col sm="2" className="d-flex justify-content-end mt-auto">
            <Button variant="primary" onClick={calculatePrevalence}>
              Calculate
            </Button>
          </Col>
        </Row>
        <div id="exposurePrevalencePlot">
          {err && (
            <div>
              <p>
                An error has occured. Check the debug section for more info.
              </p>
              <p>Error: {err}</p>
            </div>
          )}
          <div style={{ display: plotURL ? 'block' : 'none' }}>
            <Plot
              plotName={plotPath.split('/').slice(-1)[0]}
              plotURL={plotURL}
            />
          </div>
        </div>
      </Form>
      <Debug msg={debugR} />
    </div>
  );
}
