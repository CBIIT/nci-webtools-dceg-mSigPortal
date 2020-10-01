import React from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { dispatchError, dispatchExpTumor } from '../../../../services/store';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';
import Select from '../../../controls/select/select';

const { Group, Label, Control, Text } = Form;

export default function Tumor({ submitR, downloadResults }) {
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
  } = useSelector((state) => state.expTumor);
  const { displayTab, publicDataOptions } = useSelector(
    (state) => state.exploring
  );

  async function calculateR(fn, args) {
    console.log(fn);
    dispatchExpTumor({
      loading: true,
      err: false,
      debugR: '',
    });

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();

        dispatchExpTumor({
          loading: false,
          debugR: err,
        });
      } else {
        const { debugR, output } = await response.json();

        dispatchExpTumor({
          debugR: debugR,
          loading: false,
          plotPath: output.plotPath,
          txtPath: output.txtPath,
        });
        setRPlot(output.plotPath);
      }
    } catch (err) {
      dispatchError(err);
      dispatchExpTumor({ loading: false });
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
          dispatchExpTumor({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchExpTumor({ err: true, plotURL: '' });
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

    dispatchExpTumor({
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
              <Select
                id="tumorStudy"
                label="Study"
                value={study}
                options={studyOptions}
                onChange={handleStudy}
              />
            </Col>
            <Col sm="2">
              <Select
                id="tumorStrategy"
                label="Experimental Strategy"
                value={strategy}
                options={strategyOptions}
                onChange={(strategy) =>
                  dispatchExpTumor({ strategy: strategy })
                }
              />
            </Col>
            <Col sm="4">
              <Select
                id="tumorSet"
                label="Reference Signature Set"
                value={refSignatureSet}
                options={refSignatureSetOptions}
                onChange={(set) => dispatchExpTumor({ refSignatureSet: set })}
              />
            </Col>
            <Col sm="3">
              <Group controlId="tumorGenomeSize">
                <Label>Genome Size</Label>
                <Control
                  value={genomeSize}
                  onChange={(e) => {
                    dispatchExpTumor({
                      genomeSize: e.target.value,
                    });
                  }}
                />
                {/* <Text className="text-muted">(Ex. NCG>NTG)</Text> */}
              </Group>
            </Col>
            <Col sm="1" className="m-auto">
              <Button
                variant="primary"
                onClick={() => {
                  calculateR('tumorBurden', {
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
              <Plot
                plotName={plotPath.split('/').slice(-1)[0]}
                plotURL={plotURL}
                txtPath={txtPath}
              />
            </div>
          </div>
        </div>
      </Form>
      <Debug msg={debugR} />
    </div>
  );
}
