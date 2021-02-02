import React, { useEffect } from 'react';
import { Form, Row, Col, Button, Group } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { dispatchError, dispatchMsBurden } from '../../../../services/store';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';
import Select from '../../../controls/select/select';

export default function MsBurden({ calculateBurden }) {
  const { signatureNameOptions, userNameOptions, source } = useSelector(
    (state) => state.expExposure
  );
  const {
    signatureName,
    plotPath,
    plotURL,
    debugR,
    err,
    loading,
  } = useSelector((state) => state.msBurden);
  const { projectID } = useSelector((state) => state.expExposure);

  useEffect(() => {
    plotPath ? setRPlot(plotPath) : clearPlot();
  }, [plotPath, err, debugR, projectID]);

  useEffect(() => {
    if (source == 'public') {
      dispatchMsBurden({ signatureName: signatureNameOptions[0] });
    } else {
      dispatchMsBurden({ signatureName: userNameOptions[0] });
    }
  }, [signatureNameOptions, userNameOptions, source]);

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
          dispatchMsBurden({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchMsBurden({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      dispatchMsBurden({ plotURL: '' });
    }
  }

  return (
    <div>
      <Form className="p-3">
        <LoadingOverlay active={loading} />
        <Row className="">
          <Col lg="3">
            <Select
              id="acrossSignatureName"
              label="Signature Name"
              value={signatureName}
              options={
                source == 'public' ? signatureNameOptions : userNameOptions
              }
              onChange={(name) => dispatchMsBurden({ signatureName: name })}
            />
          </Col>
          <Col lg="7" />
          <Col lg="2" className="d-flex justify-content-end">
            <Button
              className="mt-auto mb-3"
              variant="primary"
              onClick={calculateBurden}
            >
              Calculate
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
              plotName="Mutational Signature Burden Across Cancer Types"
              plotURL={plotURL}
            />
          </>
        )}
        {/* <Debug msg={debugR} /> */}
      </div>
    </div>
  );
}
