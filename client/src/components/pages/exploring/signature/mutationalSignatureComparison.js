import React from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchExpMutationalSigComparison,
} from '../../../../services/store';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';
import Select from '../../../controls/select/select';

export default function MutationalSignatureProfile({ submitR }) {
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
    plotURL,
    debugR,
    err,
    loading,
  } = useSelector((state) => state.expMutationalSigComparison);
  const { displayTab, refSigData } = useSelector((state) => state.exploring);
  const { projectID } = useSelector((state) => state.visualizeResults);

  async function calculateR(fn, args) {
    dispatchExpMutationalSigComparison({
      loading: true,
      err: false,
      debugR: '',
    });

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();

        dispatchExpMutationalSigComparison({
          loading: false,
          debugR: err,
        });
      } else {
        const { debugR, output } = await response.json();
        if (Object.keys(output).length) {
          dispatchExpMutationalSigComparison({
            debugR: debugR,
            loading: false,
            plotPath: output.plotPath,
            txtPath: output.textPath,
          });
          setRPlot(output.plotPath);
        } else {
          if (plotURL) URL.revokeObjectURL(plotURL);
          dispatchExpMutationalSigComparison({
            debugR: debugR,
            loading: false,
            plotPath: '',
            plotURL: '',
            txtPath: '',
            err: true,
          });
        }
      }
    } catch (err) {
      dispatchError(err);
      dispatchExpMutationalSigComparison({ loading: false });
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
          dispatchExpMutationalSigComparison({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchExpMutationalSigComparison({ err: true, plotURL: '' });
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

    dispatchExpMutationalSigComparison({
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

    dispatchExpMutationalSigComparison({
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

    dispatchExpMutationalSigComparison({
      refSignatureSet2: set,
      signatureName2: signatureNameOptions[0],
      signatureNameOptions2: signatureNameOptions2,
    });
  }

  return (
    <div>
      <Form>
        <LoadingOverlay active={loading} />
        <div>
          <Row className="justify-content-center">
            <Col sm="3">
              <Select
                className="mb-0"
                id="mscProfileName"
                label="Profile Name"
                value={profileName}
                options={profileNameOptions}
                onChange={handleProfile}
              />
            </Col>
            <Col sm="4">
              <Select
                className="mb-0"
                id="mscRefSet1"
                label="Reference Signature Set 1"
                value={refSignatureSet1}
                options={refSignatureSetOptions1}
                onChange={handleSet1}
              />
            </Col>
            <Col sm="4">
              <Select
                className="mb-0"
                id="mscSigName1"
                label="Signature Name 1"
                value={signatureName1}
                options={signatureNameOptions1}
                onChange={(name) =>
                  dispatchExpMutationalSigComparison({
                    signatureName1: name,
                  })
                }
              />
            </Col>
            <Col sm="1" />
          </Row>
          <Row className="mt-3">
            <Col sm="3" />
            <Col sm="4">
              <Select
                className="mb-0"
                id="mscSigSet2"
                label="Signature Set 2"
                value={refSignatureSet2}
                options={refSignatureSetOptions2}
                onChange={handleSet2}
              />
            </Col>
            <Col sm="4">
              <Select
                className="mb-0"
                id="mscSetName2"
                label="Signature Name 2"
                value={signatureName2}
                options={signatureNameOptions2}
                onChange={(name) =>
                  dispatchExpMutationalSigComparison({
                    signatureName2: name,
                  })
                }
              />
            </Col>
            <Col sm="1" className="d-flex justify-content-end mt-auto">
              <Button
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
          <div id="withinPlot">
            <div style={{ display: err ? 'block' : 'none' }}>
              <p>
                An error has occured. Check the debug section for more info.
              </p>
            </div>
            <div style={{ display: plotURL ? 'block' : 'none' }}>
              <Plot
                plotName={plotPath.split('/').slice(-1)[0]}
                plotURL={plotURL}
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
