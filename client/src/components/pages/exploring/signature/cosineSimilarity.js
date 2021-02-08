import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';
import Select from '../../../controls/select/select';
import { useSelector, useDispatch } from 'react-redux';
import { actions as exploringActions } from '../../../../services/store/exploring';
import { actions as modalActions } from '../../../../services/store/modal';

const actions = { ...exploringActions, ...modalActions };

export default function MutationalSignatureProfile({ submitR }) {
  const dispatch = useDispatch();
  const exploring = useSelector((state) => state.exploring);
  const mergeExploring = (state) =>
    dispatch(actions.mergeExploring({ exploring: state }));
  const mergeSigCosineSimilarity = (state) =>
    dispatch(actions.mergeExploring({ sigCosineSimilarity: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));
  const {
    profileName,
    profileNameOptions,
    refSignatureSet1,
    refSignatureSet2,
    refSignatureSetOptions1,
    refSignatureSetOptions2,
    plotPath,
    plotURL,
    txtPath,
    debugR,
    err,
    loading,
  } = exploring.sigCosineSimilarity;
  const { displayTab, refSigData, projectID } = exploring.exploring;
  useEffect(() => {
    plotPath ? setRPlot(plotPath) : clearPlot();
  }, [plotPath]);

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      mergeSigCosineSimilarity({ plotURL: '' });
    }
  }

  async function calculateR(fn, args) {
    try {
      mergeSigCosineSimilarity({
        loading: true,
        err: false,
        plotPath: '',
        txtPath: '',
      });

      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();

        mergeSigCosineSimilarity({
          loading: false,
          debugR: err,
        });
      } else {
        const { debugR, output, projectID: id } = await response.json();
        if (Object.keys(output).length) {
          if (!projectID) mergeExploring({ projectID: id });
          mergeSigCosineSimilarity({
            debugR: debugR,
            loading: false,
            plotPath: output.plotPath,
            txtPath: output.txtPath,
          });
        } else {
          mergeSigCosineSimilarity({
            debugR: debugR,
            loading: false,
            err: true,
            plotPath: '',
            txtPath: '',
          });
        }
      }
    } catch (err) {
      mergeError(err.message);
      mergeSigCosineSimilarity({ loading: false });
    }
  }

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
          mergeSigCosineSimilarity({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        mergeError(err.message);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      mergeSigCosineSimilarity({ err: true, plotURL: '' });
    }
  }

  function handleProfile(profile) {
    let filteredData = refSigData.filter((row) => row.Profile == profile);
    const refSignatureSetOptions = [
      ...new Set(filteredData.map((row) => row.Signature_set_name)),
    ];

    mergeSigCosineSimilarity({
      profileName: profile,
      refSignatureSet1: refSignatureSetOptions[0],
      refSignatureSet2: refSignatureSetOptions[1] || refSignatureSetOptions[0],
      refSignatureSetOptions1: refSignatureSetOptions,
      refSignatureSetOptions2: refSignatureSetOptions,
    });
  }

  return (
    <div>
      <Form className="p-3">
        <LoadingOverlay active={loading} />
        <Row className="justify-content-center">
          <Col lg="3">
            <Select
              id="csProfileName"
              label="Profile Name"
              value={profileName}
              options={profileNameOptions}
              onChange={handleProfile}
            />
          </Col>
          <Col lg="4">
            <Select
              id="csRefSet1"
              label="Reference Signature Set 1"
              value={refSignatureSet1}
              options={refSignatureSetOptions1}
              onChange={(set) =>
                mergeSigCosineSimilarity({ refSignatureSet1: set })
              }
            />
          </Col>
          <Col lg="4">
            <Select
              id="rcsRefSet2"
              label="Signature Set 2"
              value={refSignatureSet2}
              options={refSignatureSetOptions2}
              onChange={(set) =>
                mergeSigCosineSimilarity({ refSignatureSet2: set })
              }
            />
          </Col>
          <Col lg="1" className="d-flex justify-content-end">
            <Button
              className="mt-auto mb-3"
              variant="primary"
              onClick={() => {
                calculateR('cosineSimilarity', {
                  profileName: profileName,
                  refSignatureSet1: refSignatureSet1,
                  refSignatureSet2: refSignatureSet2,
                });
              }}
            >
              Calculate
            </Button>
          </Col>
        </Row>
      </Form>
      <div id="withinPlot">
        <div style={{ display: err ? 'block' : 'none' }}>
          <p>An error has occured. Please verify your input.</p>
        </div>
        <hr />
        <div style={{ display: plotURL ? 'block' : 'none' }}>
          <Plot
            className="p-3"
            plotName={plotPath.split('/').slice(-1)[0]}
            plotURL={plotURL}
            txtPath={projectID + txtPath}
            maxHeight="1000px"
          />
        </div>
      </div>
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
