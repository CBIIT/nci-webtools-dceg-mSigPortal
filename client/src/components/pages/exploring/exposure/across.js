import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { dispatchError, dispatchExpAcross } from '../../../../services/store';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';
import Select from '../../../controls/select/select';

export default function Across({ calculateAcross }) {
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
  } = useSelector((state) => state.expAcross);

  useEffect(() => {
    if (plotPath) setRPlot(plotPath);
  }, [plotPath]);

  useEffect(() => {
    if (signatureNameOptions.length)
      dispatchExpAcross({ signatureName: signatureNameOptions[0] });
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

  return (
    <div>
      <Form>
        <LoadingOverlay active={loading || mainLoading} />
        <div className="px-4">
          <Row className="justify-content-center">
            <Col sm="2">
              <Select
                id="acrossSignatureName"
                label="Signature Name"
                value={signatureName}
                options={signatureNameOptions}
                onChange={(name) =>
                  dispatchExpAcross({ signatureName: name })
                }
              />
            </Col>
            <Col sm="9" />
            <Col sm="1" className="my-auto">
              <Button variant="primary" onClick={calculateAcross}>
                Calculate
              </Button>
            </Col>
          </Row>
        </div>
        <div id="acrossPlot">
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
