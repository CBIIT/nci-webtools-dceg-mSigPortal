import React, { useEffect } from 'react';
import { Form, Row, Col, Button, Group } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { dispatchError, dispatchExpAcross } from '../../../../services/store';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';
import Select from '../../../controls/select/select';

export default function Across({ calculateAcross }) {
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
  } = useSelector((state) => state.expAcross);
  const { projectID } = useSelector((state) => state.expExposure);

  useEffect(() => {
    if (plotPath) setRPlot(plotPath);
    else clearPlot();
  }, [plotPath, projectID]);

  useEffect(() => {
    if (source == 'public') {
      dispatchExpAcross({ signatureName: signatureNameOptions[0] });
    } else {
      dispatchExpAcross({ signatureName: userNameOptions[0] });
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
          dispatchExpAcross({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchExpAcross({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      dispatchExpAcross({ plotPath: '', plotURL: '' });
    }
  }

  return (
    <div>
      <Form className="p-3">
        <LoadingOverlay active={loading} />
        <p>Mutational Signature Burden Across Cancer Types</p>
        <Row className="">
          <Col lg="3">
            <Select
              id="acrossSignatureName"
              label="Signature Name"
              value={signatureName}
              options={
                source == 'public' ? signatureNameOptions : userNameOptions
              }
              onChange={(name) => dispatchExpAcross({ signatureName: name })}
            />
          </Col>
          <Col lg="7" />
          <Col lg="2" className="d-flex justify-content-end">
            <Button
              className="mt-auto mb-3"
              variant="primary"
              onClick={calculateAcross}
            >
              Calculate
            </Button>
          </Col>
        </Row>
      </Form>
      <div id="exposureAcrossPlot">
        {err && (
          <div>
            <p>An error has occured. Please verify your input.</p>
            <p>Error: {err}</p>
          </div>
        )}
        {plotURL && (
          <>
            <hr />
            <Plot
              className="p-3"
              plotName={plotPath.split('/').slice(-1)[0]}
              plotURL={plotURL}
            />
          </>
        )}
        {/* <Debug msg={debugR} /> */}
      </div>
    </div>
  );
}
