import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { actions as explorationActions } from '../../../../services/store/exploration';
import { actions as modalActions } from '../../../../services/store/modal';
import Plot from '../../../controls/plot/plot';
import Select from '../../../controls/select/select';
import Debug from '../../../controls/debug/debug';

const actions = { ...explorationActions, ...modalActions };

export default function MsBurden({ calculateBurden }) {
  const dispatch = useDispatch();
  const exploration = useSelector((state) => state.exploration);
  const {
    signatureName,
    plotPath,
    plotURL,
    debugR,
    err,
    loading,
  } = exploration.msBurden;
  const {
    projectID,
    signatureNameOptions,
    userNameOptions,
    source,
    gettingSignatureNames,
  } = exploration.exposure;
  const mergeExploration = (state) =>
    dispatch(actions.mergeExploration({ exploration: state }));
  const mergeMsBurden = (state) =>
    dispatch(actions.mergeExploration({ msBurden: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  useEffect(() => {
    plotPath ? setRPlot(plotPath) : clearPlot();
  }, [plotPath, err, debugR]);

  useEffect(() => {
    if (source == 'public') {
      mergeMsBurden({ signatureName: signatureNameOptions[0] });
    } else {
      mergeMsBurden({ signatureName: userNameOptions[0] || '' });
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
          mergeMsBurden({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        mergeError(err.message);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      mergeMsBurden({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      mergeMsBurden({ plotURL: '' });
    }
  }

  return (
    <div>
      <div className="p-3">
        <p>MS Burden: Mutational Signature Burden Across Cancer Types</p>
        <p className="m-0">
          The bar plot below illustrates the level of mutational signature
          burden across different cancer types with regard to a specific
          signature. The signature selected is one from the Reference Signature
          Set selected in the left panel. Across the top of the plot are the
          cancer types. On the y-axis is the number of mutations per Megabase
          (log10), and the x-axis denotes sample numbers. The green number is
          the number of samples with the cancer type (across the top of the
          plot), and the blue number is the number of samples in that cancer
          type that detected the signature selected.
        </p>
      </div>
      <hr />
      <Form className="p-3">
        <LoadingOverlay active={loading || gettingSignatureNames} />
        <Row className="">
          <Col lg="3">
            <Select
              className="mb-2"
              disabled={source == 'user' && !userNameOptions.length}
              id="acrossSignatureName"
              label="Signature Name"
              value={signatureName}
              options={
                source == 'public' ? signatureNameOptions : userNameOptions
              }
              onChange={(name) => mergeMsBurden({ signatureName: name })}
            />
          </Col>
          <Col />
          <Col lg="2" className="d-flex">
            <Button
              disabled={!signatureName || (source == 'user' && !projectID)}
              className="ml-auto mb-auto"
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
              title="Mutational Signature Burden Across Cancer Types"
              downloadName={plotPath.split('/').slice(-1)[0]}
              plotURL={plotURL}
            />
          </>
        )}
        {/* <Debug msg={debugR} /> */}
      </div>
    </div>
  );
}
