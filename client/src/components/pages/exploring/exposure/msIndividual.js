import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchMsIndividual,
} from '../../../../services/store';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';
import Select from '../../../controls/select/select';

const { Group } = Form;

export default function MSIndividual({ calculateIndividual }) {
  const {
    publicSampleOptions,
    userSampleOptions,
    source,
    projectID,
  } = useSelector((state) => state.expExposure);
  const { sample, plotPath, plotURL, debugR, err, loading } = useSelector(
    (state) => state.msIndividual
  );

  useEffect(() => {
    plotPath ? setRPlot(plotPath) : clearPlot();
  }, [plotPath, err, debugR, projectID]);

  // choose inital sample
  useEffect(() => {
    if (source == 'public') {
      dispatchMsIndividual({ sample: publicSampleOptions[0] });
    } else {
      dispatchMsIndividual({ sample: userSampleOptions[0] });
    }
  }, [publicSampleOptions, userSampleOptions, source]);

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
          dispatchMsIndividual({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchMsIndividual({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      dispatchMsIndividual({ plotURL: '' });
    }
  }

  return (
    <div>
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
                onChange={(name) => dispatchMsIndividual({ sample: name })}
              />
            </Group>
          </Col>
          <Col />
          <Col lg="2" className="d-flex justify-content-end">
            <Button
              disabled={
                source == 'public'
                  ? !publicSampleOptions.length
                  : !userSampleOptions.length
              }
              className="mt-auto mb-3"
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
              plotName="Mutational Signature in Individual Sample"
              plotURL={plotURL}
            />
          </>
        )}
        {/* <Debug msg={debugR} /> */}
      </div>
    </div>
  );
}
