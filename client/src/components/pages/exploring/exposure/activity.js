import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { dispatchError, dispatchExpActivity } from '../../../../services/store';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';
import Select from '../../../controls/select/select';

export default function Activity({ calculateActivity }) {
  const rootURL = window.location.pathname;
  const { signatureNameOptions, loading: mainLoading } = useSelector(
    (state) => state.expExposure
  );
  const {
    signatureName,
    plotPath,
    plotURL,
    debugR,
    err,
    loading,
  } = useSelector((state) => state.expActivity);

  useEffect(() => {
    if (plotPath) setRPlot(plotPath);
  }, [plotPath]);

  useEffect(() => {
    if (signatureNameOptions.length)
      dispatchExpActivity({ signatureName: signatureNameOptions[0] });
  }, [signatureNameOptions]);

  async function setRPlot(plotPath) {
    if (plotPath) {
      try {
        const response = await fetch(`api/results/${plotPath}`);
        if (!response.ok) {
          // console.log(await response.json());
        } else {
          const pic = await response.blob();
          const objectURL = URL.createObjectURL(pic);

          if (plotURL) URL.revokeObjectURL(plotURL);
          dispatchExpActivity({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchExpActivity({ err: true, plotURL: '' });
    }
  }

  return (
    <div>
      <Form>
        <LoadingOverlay active={loading || mainLoading} />
        <div className="px-4">
          <Row className="justify-content-center">
            <Col sm="2">
              <Select
                id="activitySignatureName"
                label="Signature Name"
                value={signatureName}
                options={signatureNameOptions}
                onChange={(name) =>
                  dispatchExpActivity({ signatureName: name })
                }
              />
            </Col>
            <Col sm="9" />
            <Col sm="1" className="my-auto">
              <Button variant="primary" onClick={calculateActivity}>
                Calculate
              </Button>
            </Col>
          </Row>
        </div>
        <div id="activityPlot">
          {err && (
            <p>An error has occured. Check the debug section for more info.</p>
          )}
          {plotURL && (
            <Plot
              plotName={plotPath.split('/').slice(-1)[0]}
              plotURL={plotURL}
            />
          )}
          <Debug msg={debugR} />
        </div>
      </Form>
    </div>
  );
}
