import React from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import Select from 'react-select';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchExpCosineSimilarity,
} from '../../../../services/store';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';

const { Group, Label } = Form;

export default function MutationalSignatureProfile({
  submitR,
  downloadResults,
}) {
  const rootURL = window.location.pathname;
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
    displayDebug,
    loading,
  } = useSelector((state) => state.expCosineSimilarity);
  const { displayTab, refSigData } = useSelector((state) => state.exploring);

  const selectFix = {
    styles: {
      menuPortal: (base) => ({ ...base, zIndex: 9999 }),
    },
    menuPortalTarget: document.body,
    getOptionLabel: (option) => option,
    getOptionValue: (option) => option,
  };

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
          txtPath: output.txtPath,
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
            <Col sm="3">
              <Group controlId="withinProfileType">
                <Label>Profile Name</Label>
                <Select
                  options={profileNameOptions}
                  value={[profileName]}
                  onChange={(profile) => handleProfile(profile)}
                  {...selectFix}
                />
              </Group>
            </Col>
            <Col sm="4">
              <Label>Reference Signature Set 1</Label>
              <Select
                options={refSignatureSetOptions1}
                value={[refSignatureSet1]}
                onChange={(set) =>
                  dispatchExpCosineSimilarity({ refSignatureSet1: set })
                }
                {...selectFix}
              />
            </Col>
            <Col sm="4">
              <Label>Signature Set 2</Label>
              <Select
                options={refSignatureSetOptions2}
                value={[refSignatureSet2]}
                onChange={(set) =>
                  dispatchExpCosineSimilarity({ refSignatureSet2: set })
                }
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
                An error has occured. Check the debug section for more info.
              </p>
            </div>
            <div style={{ display: plotURL ? 'block' : 'none' }}>
              <Plot
                plotName={plotPath.split('/').slice(-1)[0]}
                plotURL={plotURL}
                txtPath={txtPath}
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
