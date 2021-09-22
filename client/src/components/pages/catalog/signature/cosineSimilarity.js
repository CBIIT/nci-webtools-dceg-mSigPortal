import React from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';
import Select from '../../../controls/select/select';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../../services/store/catalog';
import { actions as modalActions } from '../../../../services/store/modal';

const actions = { ...catalogActions, ...modalActions };

export default function MutationalSignatureProfile({ submitR }) {
  const dispatch = useDispatch();
  const catalog = useSelector((state) => state.catalog);
  const mergeCatalog = (state) =>
    dispatch(actions.mergeCatalog({ catalog: state }));
  const mergeSigCosineSimilarity = (state) =>
    dispatch(actions.mergeCatalog({ sigCosineSimilarity: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));
  const {
    profileName,
    profileNameOptions,
    rsSet1,
    rsSet2,
    rsSetOptions1,
    rsSetOptions2,
    plotPath,
    txtPath,
    debugR,
    err,
    loading,
  } = catalog.sigCosineSimilarity;
  const { refSigData, projectID } = catalog.catalog;

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
        const { stdout, output, projectID: id } = await response.json();
        if (output.plotPath) {
          if (!projectID) mergeCatalog({ projectID: id });
          mergeSigCosineSimilarity({
            loading: false,
            plotPath: output.plotPath,
            txtPath: output.txtPath,
          });
        } else {
          mergeSigCosineSimilarity({
            loading: false,
            err: output.error || output.uncaughtError || true,
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

  function handleProfile(profile) {
    let filteredData = refSigData.filter((row) => row.Profile == profile);
    const rsSetOptions = [
      ...new Set(filteredData.map((row) => row.Signature_set_name)),
    ];

    mergeSigCosineSimilarity({
      profileName: profile,
      rsSet1: rsSetOptions[0],
      rsSet2: rsSetOptions[1] || rsSetOptions[0],
      rsSetOptions1: rsSetOptions,
      rsSetOptions2: rsSetOptions,
    });
  }

  return (
    <div style={{ minHeight: '500px' }}>
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
              value={rsSet1}
              options={rsSetOptions1}
              onChange={(set) => mergeSigCosineSimilarity({ rsSet1: set })}
            />
          </Col>
          <Col lg="4">
            <Select
              id="rcsRefSet2"
              label="Reference Signature Set 2"
              value={rsSet2}
              options={rsSetOptions2}
              onChange={(set) => mergeSigCosineSimilarity({ rsSet2: set })}
            />
          </Col>
          <Col lg="1" className="d-flex justify-content-end">
            <Button
              className="mt-auto mb-3"
              variant="primary"
              onClick={() => {
                calculateR('cosineSimilarity', {
                  profileName: profileName,
                  rsSet1: rsSet1,
                  rsSet2: rsSet2,
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
        {plotPath && (
          <>
            <Plot
              className="p-3"
              downloadName={plotPath.split('/').slice(-1)[0]}
              plotPath={`api/results/${plotPath}`}
              txtPath={projectID + txtPath}
              height="1000px"
              title="Cosine Similarity Among Mutational Signatures Between Reference Signature Sets"
            />
            <p className="p-3">
              The Cosine Similarity Among Mutational Signatures Between
              Reference Signature Sets plot highlights cosine similarity between
              two mutational signature sets given a profile type. Along the
              bottom of the heatmap are the signatures within the reference
              signature set selected. Along the side of the heatmap are the
              signatures within the second reference signature selected.
            </p>
          </>
        )}
      </div>
    </div>
  );
}
