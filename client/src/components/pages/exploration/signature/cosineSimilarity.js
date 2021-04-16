import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';
import Select from '../../../controls/select/select';
import { useSelector, useDispatch } from 'react-redux';
import { actions as explorationActions } from '../../../../services/store/exploration';
import { actions as modalActions } from '../../../../services/store/modal';

const actions = { ...explorationActions, ...modalActions };

export default function MutationalSignatureProfile({ submitR }) {
  const dispatch = useDispatch();
  const exploration = useSelector((state) => state.exploration);
  const mergeExploration = (state) =>
    dispatch(actions.mergeExploration({ exploration: state }));
  const mergeSigCosineSimilarity = (state) =>
    dispatch(actions.mergeExploration({ sigCosineSimilarity: state }));
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
  } = exploration.sigCosineSimilarity;
  const { displayTab, refSigData, projectID } = exploration.exploration;
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
          if (!projectID) mergeExploration({ projectID: id });
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
      <div className="p-3">
        <p>
          Cosine similarity, which is a measure of the similarity of signatures,
          can be helpful to compare one signature to another. In the case below,
          you are able to compare two reference mutational signatures. Simply
          use the dropdown menus to enter a “Profile Name”, a “Reference
          Signature Set 1”, and “Reference Signature Set 2”. This will compare
          the mutational profile entered between the two reference signature
          sets. Click here to learn more about cosine.
        </p>
      </div>
      <hr />
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
              label="Reference Signature Set 2"
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
        <div style={{ display: err ? 'block' : 'none' }} className="p-3">
          <p>An error has occured. Please verify your input.</p>
        </div>
        <hr />
        <div style={{ display: plotURL ? 'block' : 'none' }}>
          <Plot
            className="p-3"
            downloadName={plotPath.split('/').slice(-1)[0]}
            plotURL={plotURL}
            txtPath={projectID + txtPath}
            maxHeight="1000px"
            title="Cosine Similarity Among Mutational Signatures Between Reference Signature Sets"
          />
          <p className="p-3">
            The Cosine Similarity Among Mutational Signatures Between Reference
            Signature Sets plot highlights cosine similarity between two
            mutational signature sets given a profile type. Along the bottom of
            the heatmap are the signatures within the reference signature set
            selected. Along the side of the heatmap are the signatures within
            the second reference signature selected.
          </p>
        </div>
      </div>
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
