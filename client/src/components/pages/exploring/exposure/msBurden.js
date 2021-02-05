import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { actions as exploringActions } from '../../../../services/store/exploring';
import { actions as modalActions } from '../../../../services/store/modal';
import Plot from '../../../controls/plot/plot';
import Select from '../../../controls/select/select';
import Debug from '../../../controls/debug/debug';

const actions = { ...exploringActions, ...modalActions };

export default function MsBurden({ calculateBurden }) {
  const dispatch = useDispatch();
  const exploring = useSelector((state) => state.exploring);
  const {
    signatureName,
    plotPath,
    plotURL,
    debugR,
    err,
    loading,
  } = exploring.msBurden;
  const {
    projectID,
    signatureNameOptions,
    userNameOptions,
    source,
  } = exploring.exposure;
  const mergeExploring = (state) =>
    dispatch(actions.mergeExploring({ exploring: state }));
  const mergeMsBurden = (state) =>
    dispatch(actions.mergeExploring({ msBurden: state }));
  const mergeError = (state) => dispatch(actions.mergeModal({ error: state }));

  useEffect(() => {
    plotPath ? setRPlot(plotPath) : clearPlot();
  }, [plotPath, err, debugR, projectID]);

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
        mergeError({ visible: true, message: err.message });
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
      <Form className="p-3">
        <LoadingOverlay active={loading} />
        <Row className="">
          <Col lg="3">
            <Select
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
          <Col lg="7" />
          <Col lg="2" className="d-flex justify-content-end">
            <Button
              disabled={!signatureName}
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
