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
const { Group } = Form;

export default function MSIndividual({ calculateIndividual }) {
  const dispatch = useDispatch();
  const exploration = useSelector((state) => state.exploration);
  const {
    sample,
    plotPath,
    plotURL,
    debugR,
    err,
    loading,
  } = exploration.msIndividual;
  const {
    projectID,
    publicSampleOptions,
    userSampleOptions,
    source,
  } = exploration.exposure;
  const mergeExploration = (state) =>
    dispatch(actions.mergeExploration({ exploration: state }));
  const mergeMsIndividual = (state) =>
    dispatch(actions.mergeExploration({ msIndividual: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  useEffect(() => {
    plotPath ? setRPlot(plotPath) : clearPlot();
  }, [plotPath, err, debugR]);

  // choose inital sample
  useEffect(() => {
    if (source == 'public') {
      mergeMsIndividual({ sample: publicSampleOptions[0] });
    } else {
      mergeMsIndividual({ sample: userSampleOptions[0] || '' });
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
          mergeMsIndividual({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        mergeError(err.message);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      mergeMsIndividual({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      mergeMsIndividual({ plotURL: '' });
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
                onChange={(name) => mergeMsIndividual({ sample: name })}
              />
            </Group>
          </Col>
          <Col />
          <Col lg="2" className="d-flex">
            <Button
              className="ml-auto mb-auto"
              disabled={!sample || (source == 'user' && !projectID)}
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
              title="Mutational Signature in Individual Sample"
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
