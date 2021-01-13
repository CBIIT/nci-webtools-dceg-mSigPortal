import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchExpCosineSimilarity,
  dispatchExploring,
} from '../../../../services/store';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';
import Select from '../../../controls/select/select';

export default function MutationalSignatureProfile({ submitR }) {
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
  } = useSelector((state) => state.expCosineSimilarity);
  const { displayTab, refSigData, projectID } = useSelector(
    (state) => state.exploring
  );

  useEffect(() => {
    if (plotPath) setRPlot(plotPath);
    else clearPlot();
  }, [plotPath]);

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      dispatchExpCosineSimilarity({ plotPath: '', plotURL: '' });
    }
  }

  async function calculateR(fn, args) {
    dispatchExpCosineSimilarity({
      loading: true,
      err: false,
      debugR: '',
    });

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();

        dispatchExpCosineSimilarity({
          loading: false,
          debugR: err,
        });
      } else {
        const { debugR, output, projectID: id } = await response.json();
        if (Object.keys(output).length) {
          if (!projectID) dispatchExploring({ projectID: id });
          dispatchExpCosineSimilarity({
            debugR: debugR,
            loading: false,
            plotPath: output.plotPath,
            txtPath: output.txtPath,
          });
        } else {
          dispatchExpCosineSimilarity({
            debugR: debugR,
            loading: false,
            err: true,
            plotPath: '',
            txtPath: '',
          });
        }
      }
    } catch (err) {
      dispatchError(err);
      dispatchExpCosineSimilarity({ loading: false });
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
          dispatchExpCosineSimilarity({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchExpCosineSimilarity({ err: true, plotURL: '' });
    }
  }

  function handleProfile(profile) {
    let filteredData = refSigData.filter((row) => row.Profile == profile);
    const refSignatureSetOptions = [
      ...new Set(filteredData.map((row) => row.Signature_set_name)),
    ];

    dispatchExpCosineSimilarity({
      profileName: profile,
      refSignatureSet1: refSignatureSetOptions[0],
      refSignatureSet2: refSignatureSetOptions[1] || refSignatureSetOptions[0],
      refSignatureSetOptions1: refSignatureSetOptions,
      refSignatureSetOptions2: refSignatureSetOptions,
    });
  }

  return (
    <div>
      <Form>
        <LoadingOverlay active={loading} />
        <div>
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
                  dispatchExpCosineSimilarity({ refSignatureSet1: set })
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
                  dispatchExpCosineSimilarity({ refSignatureSet2: set })
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
          <div id="withinPlot">
            <div style={{ display: err ? 'block' : 'none' }}>
              <p>
                An error has occured. Please verify your input.
              </p>
            </div>
            <div style={{ display: plotURL ? 'block' : 'none' }}>
              <Plot
                plotName={plotPath.split('/').slice(-1)[0]}
                plotURL={plotURL}
                txtPath={projectID + txtPath}
                maxHeight="1000px"
              />
            </div>
          </div>
        </div>
      </Form>
      <Debug msg={debugR} />
    </div>
  );
}
