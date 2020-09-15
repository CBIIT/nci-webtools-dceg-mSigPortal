import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import Select from 'react-select';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchExpCosineSimilarity,
} from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

const { Group, Label } = Form;

export default function MutationalSignatureProfile({
  submitR,
  downloadResults,
  getRefSigOptions,
}) {
  const rootURL = window.location.pathname;
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
    txtPath,
    debugR,
    err,
    displayDebug,
    loading,
  } = useSelector((state) => state.expCosineSimilarity);

  const selectFix = {
    styles: {
      menuPortal: (base) => ({ ...base, zIndex: 9999 }),
    },
    menuPortalTarget: document.body,
    getOptionLabel: (option) => option,
    getOptionValue: (option) => option,
  };

  const { displayTab } = useSelector((state) => state.exploring);

  useEffect(() => {
    if (!loading && !plotPath && displayTab == 'cosineSimilarity') {
      calculateR('cosineSimilarity', {});
    }
  }, [plotPath, displayTab]);

  async function calculateR(fn, args) {
    console.log(fn);
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
        const { debugR, output } = await response.json();

        dispatchExpCosineSimilarity({
          debugR: debugR,
          loading: false,
          plotPath: output.plotPath,
          txtPath: output.textPath,
        });
        setRPlot(output.plotPath);
      }
    } catch (err) {
      dispatchError(err);
      dispatchExpCosineSimilarity({ loading: false });
    }
  }

  async function setRPlot(plotPath) {
    if (plotPath) {
      try {
        const response = await fetch(`${rootURL}getSVG`, {
          method: 'POST',
          headers: {
            Accept: 'image/svg',
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({ path: plotPath }),
        });
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

  // get Signature Reference Sets for dropdown options
  async function getSignatureSet(profileType) {
    if (profileType) {
      dispatchExpCosineSimilarity({ loading: true });
      try {
        const response = await getRefSigOptions(profileType);

        if (response.ok) {
          const refSignatureSetOptions = await response.json();

          dispatchExpCosineSimilarity({
            refSignatureSetOptions: refSignatureSetOptions,
            refSignatureSet: refSignatureSetOptions[0],
            loading: false,
          });
        } else {
          dispatchError(await response.json());
          dispatchExpCosineSimilarity({ loading: false });
        }
      } catch (err) {
        dispatchError(err);
        dispatchExpCosineSimilarity({ loading: false });
      }
    }
  }

  function handleProfileType(profileType) {
    const matrixList = [];
    const matrixOptions = [
      ...new Set(
        matrixList
          .filter((matrix) => matrix.Profile_Type == profileType)
          .map((matrix) => matrix.Matrix_Size)
      ),
    ];

    dispatchExpCosineSimilarity({
      profile: profileType,
      matrix: matrixOptions[0],
      matrixOptions: matrixOptions,
    });
  }

  return (
    <div>
      <Form>
        <LoadingOverlay active={loading} />
        <div>
          <Row className="justify-content-center">
            <Col sm="3">
              <Group controlId="withinProfileType">
                <Label>Profile Name</Label>
                <Select
                  options={profileNameOptions}
                  value={[profileName]}
                  // onChange={(profile) => handleProfileName(profile)}
                  {...selectFix}
                />
              </Group>
            </Col>
            <Col sm="4">
              <Label>Reference Signature Set 1</Label>
              <Select
                options={refSignatureSetOptions1}
                value={[refSignatureSet1]}
                // onChange={(set) =>
                //   dispatchExpCosineSimilarity({
                //     refSignatureSet1: set,
                //   })
                // }
                {...selectFix}
              />
            </Col>
            <Col sm="4">
              <Label>Signature Name 1</Label>
              <Select
                options={signatureNameOptions1}
                value={[signatureName1]}
                // onChange={(name) =>
                //   dispatchExpCosineSimilarity({
                //     signatureName1: name,
                //   })
                // }
                {...selectFix}
              />
            </Col>
            <Col sm="1" />
          </Row>
          <Row>
            <Col sm="3" />
            <Col sm="4">
              <Label>Signature Set 2</Label>
              <Select
                options={refSignatureSet2}
                value={[refSignatureSetOptions2]}
                onChange={(set) =>
                  dispatchExpCosineSimilarity({
                    refSignatureSet2: set,
                  })
                }
                {...selectFix}
              />
            </Col>
            <Col sm="4">
              <Label>Signature Name 2</Label>
              <Select
                options={signatureNameOptions2}
                value={[signatureName2]}
                // onChange={(name) =>
                //   dispatchExpCosineSimilarity({
                //     signatureName1: name,
                //   })
                // }
                {...selectFix}
              />
            </Col>
            <Col sm="1" className="m-auto">
              <Button
                variant="primary"
                onClick={() => {
                  calculateR('cosineSimilarity', {
                    profileName: profileName,
                    refSignatureSet1: refSignatureSet1,
                    signatureName1: signatureNameOptions1,
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
              <div className="d-flex">
                <a
                  className="px-2 py-1"
                  href={plotURL}
                  download={plotURL.split('/').slice(-1)[0]}
                >
                  Download Plot
                </a>
                <span className="ml-auto">
                  <Button
                    className="px-2 py-1"
                    variant="link"
                    onClick={() => downloadResults(txtPath)}
                  >
                    Download Results
                  </Button>
                </span>
              </div>
              <div className="p-2 border rounded">
                <Row>
                  <Col>
                    <img className="w-100 my-4 h-1000" src={plotURL}></img>
                  </Col>
                </Row>
              </div>
            </div>
          </div>
        </div>
      </Form>
      <Button
        variant="link"
        className="p-0 mt-5"
        onClick={() =>
          dispatchExpCosineSimilarity({
            displayDebug: !displayDebug,
          })
        }
      >
        R Debug
      </Button>
      <pre
        className="border rounded p-1 "
        style={{ display: displayDebug ? 'block' : 'none' }}
      >
        <div className="border">
          {Array.isArray(debugR) ? (
            debugR.map((line, index) => {
              return (
                <p key={index} className="m-0">
                  [{index}] {line}
                </p>
              );
            })
          ) : (
            <p>{debugR}</p>
          )}
        </div>
      </pre>
    </div>
  );
}
