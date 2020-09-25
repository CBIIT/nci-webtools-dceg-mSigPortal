import React from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import Select from 'react-select';
import { useSelector } from 'react-redux';
import { dispatchError, dispatchExpActivity } from '../../../../services/store';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

const { Group, Label, Control, Text } = Form;

export default function Activity({ submitR, downloadResults }) {
  const rootURL = window.location.pathname;
  const {
    study,
    studyOptions,
    strategy,
    strategyOptions,
    refSignatureSet,
    refSignatureSetOptions,
    genomeSize,
    plotPath,
    plotURL,
    txtPath,
    debugR,
    err,
    displayDebug,
    loading,
  } = useSelector((state) => state.expActivity);
  const { displayTab, publicDataOptions } = useSelector(
    (state) => state.exploring
  );

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
    dispatchExpActivity({
      loading: true,
      err: false,
      debugR: '',
    });

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();

        dispatchExpActivity({
          loading: false,
          debugR: err,
        });
      } else {
        const { debugR, output } = await response.json();

        dispatchExpActivity({
          debugR: debugR,
          loading: false,
          plotPath: output.plotPath,
          txtPath: output.txtPath,
        });
        setRPlot(output.plotPath);
      }
    } catch (err) {
      dispatchError(err);
      dispatchExpActivity({ loading: false });
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
          dispatchExpActivity({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchExpActivity({ err: true, plotURL: '' });
    }
  }

  function handleStudy(study) {
    const strategyOptions = [
      ...new Set(
        publicDataOptions
          .filter((data) => data.Study == study)
          .map((data) => data.Dataset)
      ),
    ];

    dispatchExpActivity({
      study: study,
      strategy: strategyOptions[0],
      strategyOptions: strategyOptions,
    });
  }

  return (
    <div>
      <Form>
        <LoadingOverlay active={loading} />
        <div>
          <Row className="justify-content-center">
            <Col sm="2">
              <Group controlId="withinProfileType">
                <Label>Study</Label>
                <Select
                  options={studyOptions}
                  value={[study]}
                  onChange={(study) => handleStudy(study)}
                  {...selectFix}
                />
              </Group>
            </Col>
            <Col sm="2">
              <Label>Experimental Strategy</Label>
              <Select
                options={strategyOptions}
                value={[strategy]}
                onChange={(strategy) =>
                  dispatchExpActivity({ strategy: strategy })
                }
                {...selectFix}
              />
            </Col>
            <Col sm="4">
              <Label>Reference Signature Set</Label>
              <Select
                options={refSignatureSetOptions}
                value={[refSignatureSet]}
                onChange={(set) =>
                  dispatchExpActivity({ refSignatureSet: set })
                }
                {...selectFix}
              />
            </Col>
            <Col sm="3">
              <Label>Genome Size</Label>
              <Control
                value={genomeSize}
                onChange={(e) => {
                  dispatchExpActivity({
                    genomeSize: e.target.value,
                  });
                }}
              />
              {/* <Text className="text-muted">(Ex. NCG>NTG)</Text> */}
            </Col>
            <Col sm="1" className="m-auto">
              <Button
                variant="primary"
                onClick={() => {
                  calculateR('cosineSimilarity', {
                    study: study,
                    strategy: strategy,
                    refSignatureSet: refSignatureSet,
                    genomeSize: parseFloat(genomeSize),
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
                  download={plotPath.split('/').slice(-1)[0]}
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
          dispatchExpActivity({
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
