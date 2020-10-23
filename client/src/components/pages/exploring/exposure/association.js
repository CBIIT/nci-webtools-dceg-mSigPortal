import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchExpAssociation,
} from '../../../../services/store';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';
import Select from '../../../controls/select/select';

export default function Association({ submitR }) {
  const rootURL = window.location.pathname;
  const { signatureNameOptions, loading: mainLoading } = useSelector(
    (state) => state.expExposure
  );
  const {
    signatureName1,
    signatureName2,
    plotPath,
    plotURL,
    debugR,
    err,
    loading,
  } = useSelector((state) => state.expAssociation);

  useEffect(() => {
    if (plotPath) setRPlot(plotPath);
  }, [plotPath]);

  useEffect(() => {
    if (signatureNameOptions.length)
      dispatchExpAssociation({
        signatureName1: signatureNameOptions[0],
        signatureName2: signatureNameOptions[1] || signatureNameOptions[0],
      });
  }, [signatureNameOptions]);

  // async function calculateR(fn, args) {
  //   console.log(fn);
  //   dispatchExpAssociation({
  //     loading: true,
  //     err: false,
  //     debugR: '',
  //   });

  //   try {
  //     const response = await submitR(fn, args);
  //     if (!response.ok) {
  //       const err = await response.json();

  //       dispatchExpAssociation({
  //         loading: false,
  //         debugR: err,
  //       });
  //     } else {
  //       const { debugR, output } = await response.json();

  //       dispatchExpAssociation({
  //         debugR: debugR,
  //         loading: false,
  //         plotPath: output.plotPath,
  //         txtPath: output.txtPath,
  //       });
  //       setRPlot(output.plotPath);
  //     }
  //   } catch (err) {
  //     dispatchError(err);
  //     dispatchExpAssociation({ loading: false });
  //   }
  // }

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
          dispatchExpAssociation({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchExpAssociation({ err: true, plotURL: '' });
    }
  }

  return (
    <div>
      <Form>
        <LoadingOverlay active={loading || mainLoading} />
        <div>
          <Row className="justify-content-center">
            <Col sm="2">
              <Select
                id="associationSignatureName1"
                label="Signature Name 1"
                value={signatureName1}
                options={signatureNameOptions}
                onChange={(name) =>
                  dispatchExpAssociation({ signatureName1: name })
                }
              />
            </Col>
            <Col sm="2">
              <Select
                id="associationSignatureName2"
                label="Signature Name 2"
                value={signatureName2}
                options={signatureNameOptions}
                onChange={(name) =>
                  dispatchExpAssociation({ signatureName2: name })
                }
              />
            </Col>
            <Col sm="7" />
            <Col sm="1" className="my-auto">
              <Button
                variant="primary"
                onClick={() => {
                  // calculateR('cosineSimilarity', {
                  //   study: study,
                  //   strategy: strategy,
                  //   refSignatureSet: refSignatureSet,
                  //   genomeSize: parseFloat(genomeSize),
                  // });
                }}
              >
                Calculate
              </Button>
            </Col>
          </Row>
          <div id="associationPlot">
            {err && (
              <p>
                An error has occured. Check the debug section for more info.
              </p>
            )}
            {plotURL && (
              <Plot
                plotName={plotPath.split('/').slice(-1)[0]}
                plotURL={plotURL}
              />
            )}
            <Debug msg={debugR} />
          </div>
        </div>
      </Form>
    </div>
  );
}
