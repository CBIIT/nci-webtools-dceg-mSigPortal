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
  const mergeSigMutationalSigComparison = (state) =>
    dispatch(actions.mergeExploration({ sigMutationalSigComparison: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const {
    profileName,
    profileNameOptions,
    signatureName1,
    signatureNameOptions1,
    signatureName2,
    signatureNameOptions2,
    refSignatureSet1,
    refSignatureSetOptions1,
    refSignatureSet2,
    refSignatureSetOptions2,
    plotPath,
    debugR,
    err,
    loading,
  } = exploration.sigMutationalSigComparison;
  const { displayTab, refSigData, projectID } = exploration.exploration;

  async function calculateR(fn, args) {
    mergeSigMutationalSigComparison({
      loading: true,
      err: false,
      plotPath: '',
      txtPath: '',
    });

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();

        mergeSigMutationalSigComparison({
          loading: false,
          debugR: err,
        });
      } else {
        const { debugR, output, projectID: id } = await response.json();
        if (Object.keys(output).length) {
          if (!projectID) mergeExploration({ projectID: id });

          mergeSigMutationalSigComparison({
            debugR: debugR,
            loading: false,
            plotPath: output.plotPath,
            txtPath: output.textPath,
          });
        } else {
          mergeSigMutationalSigComparison({
            debugR: debugR,
            loading: false,
            plotPath: '',
            txtPath: '',
            err: true,
          });
        }
      }
    } catch (err) {
      mergeError(err.message);
      mergeSigMutationalSigComparison({ loading: false });
    }
  }

  function handleProfile(profile) {
    let filteredData = refSigData.filter((row) => row.Profile == profile);
    const refSignatureSetOptions = [
      ...new Set(filteredData.map((row) => row.Signature_set_name)),
    ];
    const refSignatureSet1 = refSignatureSetOptions[0];
    const refSignatureSet2 =
      refSignatureSetOptions[1] || refSignatureSetOptions[0];
    const signatureNameOptions1 = [
      ...new Set(
        filteredData
          .filter((row) => row.Signature_set_name == refSignatureSet1)
          .map((row) => row.Signature_name)
      ),
    ];
    const signatureNameOptions2 = [
      ...new Set(
        filteredData
          .filter((row) => row.Signature_set_name == refSignatureSet2)
          .map((row) => row.Signature_name)
      ),
    ];

    mergeSigMutationalSigComparison({
      profileName: profile,
      refSignatureSet1: refSignatureSet1,
      refSignatureSet2: refSignatureSet2,
      refSignatureSetOptions1: refSignatureSetOptions,
      refSignatureSetOptions2: refSignatureSetOptions,
      signatureName1: signatureNameOptions1[0],
      signatureName2: signatureNameOptions2[0],
      signatureNameOptions1: signatureNameOptions1,
      signatureNameOptions2: signatureNameOptions2,
    });
  }

  function handleSet1(set) {
    let filteredData = refSigData.filter(
      (row) => row.Profile == profileName && row.Signature_set_name == set
    );

    const signatureNameOptions1 = [
      ...new Set(filteredData.map((row) => row.Signature_name)),
    ];

    mergeSigMutationalSigComparison({
      refSignatureSet1: set,
      signatureName1: signatureNameOptions1[0],
      signatureNameOptions1: signatureNameOptions1,
    });
  }

  function handleSet2(set) {
    let filteredData = refSigData.filter(
      (row) => row.Profile == profileName && row.Signature_set_name == set
    );

    const signatureNameOptions = [
      ...new Set(filteredData.map((row) => row.Signature_name)),
    ];

    mergeSigMutationalSigComparison({
      refSignatureSet2: set,
      signatureName2: signatureNameOptions[0],
      signatureNameOptions2: signatureNameOptions,
    });
  }

  return (
    <div>
      <p className="p-3">
        Below you can observe mutational signature comparisons signatures found
        in two different reference signature sets. Use the dropdown menus to
        input a “Profile Name”, two Reference Signature Sets, and two Signatures
        within the selected Signature Sets. This will allow for a calculation as
        to the mutational profile of each signature given the profile type, as
        well as the difference between the two mutational profiles.
      </p>
      <hr />
      <Form className="p-3">
        <LoadingOverlay active={loading} />
        <Row className="justify-content-center">
          <Col lg="3">
            <Select
              id="mscProfileName"
              label="Profile Name"
              value={profileName}
              options={profileNameOptions}
              onChange={handleProfile}
            />
          </Col>
          <Col lg="4">
            <Select
              id="mscRefSet1"
              label="Reference Signature Set 1"
              value={refSignatureSet1}
              options={refSignatureSetOptions1}
              onChange={handleSet1}
            />
          </Col>
          <Col lg="4">
            <Select
              id="mscSigName1"
              label="Signature Name 1"
              value={signatureName1}
              options={signatureNameOptions1}
              onChange={(name) =>
                mergeSigMutationalSigComparison({
                  signatureName1: name,
                })
              }
            />
          </Col>
          <Col lg="1" />
        </Row>
        <Row>
          <Col lg="3" />
          <Col lg="4">
            <Select
              id="mscSigSet2"
              label="Reference Signature Set 2"
              value={refSignatureSet2}
              options={refSignatureSetOptions2}
              onChange={handleSet2}
            />
          </Col>
          <Col lg="4">
            <Select
              id="mscSetName2"
              label="Signature Name 2"
              value={signatureName2}
              options={signatureNameOptions2}
              onChange={(name) =>
                mergeSigMutationalSigComparison({
                  signatureName2: name,
                })
              }
            />
          </Col>
          <Col lg="1" className="d-flex justify-content-end">
            <Button
              className="mt-auto mb-3"
              variant="primary"
              onClick={() => {
                calculateR('mutationalSignatureComparison', {
                  profileName: profileName,
                  refSignatureSet1: refSignatureSet1,
                  signatureName1: signatureName1,
                  refSignatureSet2: refSignatureSet2,
                  signatureName2: signatureName2,
                });
              }}
            >
              Calculate
            </Button>
          </Col>
        </Row>
      </Form>

      <div id="mutationalSignatureComparison">
        <div style={{ display: err ? 'block' : 'none' }} className="p-3">
          <p>An error has occured. Please verify your input.</p>
        </div>
        {plotPath && (
          <>
            <hr />
            <Plot
              className="p-3"
              downloadName={plotPath.split('/').slice(-1)[0]}
              plotPath={`api/results/${projectID}${plotPath}`}
              maxHeight="1000px"
            />
            <div className="p-3">
              <p>
                The plot generated shows the mutational profile for the
                “Signature Name 1” selected, the mutational profile for the
                “Signature Name 2” selected, and the difference between them.
                Also at the top of the plot are measurements for RSS and cosine
                similarity.
              </p>
              <p>
                RSS is the Residual Sum of Squares. It measures the discrepancy
                between two profiles. Cosine similarity is how similar the
                mutational profiles are to one another. For additional
                information about RSS and cosine similarity, click here.
              </p>
            </div>
          </>
        )}
      </div>
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
