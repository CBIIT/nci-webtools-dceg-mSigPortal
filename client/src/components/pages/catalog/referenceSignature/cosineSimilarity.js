import React from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import SvgContainer from '../../../controls/svgContainer/svgContainer';
import CustomSelect from '../../../controls/select/select-old';
import Description from '../../../controls/description/description';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../../services/store/catalog';
import { actions as modalActions } from '../../../../services/store/modal';
import { NavHashLink } from 'react-router-hash-link';

const actions = { ...catalogActions, ...modalActions };

export default function MutationalSignatureProfile({ submitR }) {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.catalog);
  const mergeRS = (state) =>
    dispatch(actions.mergeCatalog({ referenceSignature: state }));
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
  } = store.sigCosineSimilarity;
  const { refSigData, projectID } = store.referenceSignature;

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
          if (!projectID) mergeRS({ projectID: id });
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
    <div>
      <Description
        className="p-3"
        less="Cosine similarity is a measure of the similarity of two signature matrices, which can be helpful to compare two mutational profiles or signatures. Below you can explore cosine similarity between two reference mutational signature sets."
        more={
          <span>
            Use the dropdown menus to enter a [Profile Name], [Reference
            Signature Set 1], and [Reference Signature Set 2]. Click{' '}
            <NavHashLink to="/faq#cosine-similarity">here</NavHashLink> to learn
            more about cosine similarity.
          </span>
        }
      />
      <hr />
      <Form className="p-3">
        <LoadingOverlay active={loading} />
        <Row className="">
          <Col lg="auto">
            <CustomSelect
              id="csProfileName"
              label="Profile Name"
              value={profileName}
              options={profileNameOptions}
              onChange={handleProfile}
            />
          </Col>
          <Col lg="auto">
            <CustomSelect
              id="csRefSet1"
              label="Reference Signature Set 1"
              value={rsSet1}
              options={rsSetOptions1}
              onChange={(set) => mergeSigCosineSimilarity({ rsSet1: set })}
            />
          </Col>
          <Col lg="auto">
            <CustomSelect
              id="rcsRefSet2"
              label="Reference Signature Set 2"
              value={rsSet2}
              options={rsSetOptions2}
              onChange={(set) => mergeSigCosineSimilarity({ rsSet2: set })}
            />
          </Col>
          <Col lg="auto" className="d-flex justify-content-end">
            <Button
              className="mt-auto mb-3"
              variant="primary"
              onClick={() => {
                calculateR('rsCosineSimilarity', {
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
            <SvgContainer
              className="p-3"
              downloadName={plotPath.split('/').slice(-1)[0]}
              plotPath={`web/results/${plotPath}`}
              txtPath={`web/results/${txtPath}`}
              height="1000px"
              title="Cosine Similarity Among Mutational Signatures Between Two Reference Signature Sets"
            />
            <p className="p-3">
              The heatmap above shows the cosine similarities between two
              mutational signature sets given a profile type. The text on the
              bottom and left show the signature names of the two selected
              reference signature sets.
            </p>
          </>
        )}
      </div>
    </div>
  );
}
