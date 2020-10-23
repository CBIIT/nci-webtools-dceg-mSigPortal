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

const { Group, Label, Control, Text } = Form;

export default function Tumor({ submitR }) {
  const rootURL = window.location.pathname;
  const { loading: mainLoading } = useSelector((state) => state.expExposure);
  const { mutation, plotPath, plotURL, debugR, err, loading } = useSelector(
    (state) => state.expPrevalence
  );

  useEffect(() => {
    if (plotPath) setRPlot(plotPath);
  }, [plotPath]);

  async function calculateR(fn, args) {
    console.log(fn);
    dispatchExpPrevalence({
      loading: true,
      err: false,
      debugR: '',
    });

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();

        dispatchExpPrevalence({
          loading: false,
          debugR: err,
        });
      } else {
        const { debugR, output } = await response.json();

        dispatchExpPrevalence({
          debugR: debugR,
          loading: false,
          plotPath: output.plotPath,
          txtPath: output.txtPath,
        });
        setRPlot(output.plotPath);
      }
    } catch (err) {
      dispatchError(err);
      dispatchExpPrevalence({ loading: false });
    }
  }

  async function setRPlot(plotPath) {
    if (plotPath) {
      try {
        const response = await fetch(`${rootURL}results/${plotPath}`);
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

  return (
    <div>
      <Form>
        <LoadingOverlay active={loading || mainLoading} />
        <div>
          <Row className="justify-content-center">
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
            <Col sm="7" />
            <Col sm="1" className="m-auto">
              <Button variant="primary" onClick={() => {}}>
                Calculate
              </Button>
            </Col>
          </Row>
          <div id="withinPlot">
            <div style={{ display: err ? 'block' : 'none' }}>
              <p>
                An error has occured. Check the debug section for more info.
              </p>
            </div>
            <div style={{ display: plotURL ? 'block' : 'none' }}>
              <Plot
                plotName={plotPath.split('/').slice(-1)[0]}
                plotURL={plotURL}
              />
            </div>
          </div>
        </div>
      </Form>
      <Debug msg={debugR} />
    </div>
  );
}
