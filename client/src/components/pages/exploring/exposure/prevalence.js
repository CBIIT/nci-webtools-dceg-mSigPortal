import React from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchExpPrevalence,
} from '../../../../services/store';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';
import Select from '../../../controls/select/select';

const { Group, Label, Control, Text } = Form;

export default function Tumor({ submitR }) {
  const rootURL = window.location.pathname;
  const {
    study,
    studyOptions,
    strategy,
    strategyOptions,
    refSignatureSet,
    refSignatureSetOptions,
    genomeSize,
    mutation,
    plotPath,
    plotURL,
    txtPath,
    debugR,
    err,
    displayDebug,
    loading,
  } = useSelector((state) => state.expPrevalence);
  const { displayTab, publicDataOptions } = useSelector(
    (state) => state.exploring
  );

  async function calculateR(fn, args) {
    console.log(fn);
    dispatchExpPrevalence({
      loading: true,
      err: false,
      debugR: '',
    });

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();

        dispatchExpPrevalence({
          loading: false,
          debugR: err,
        });
      } else {
        const { debugR, output } = await response.json();

        dispatchExpPrevalence({
          debugR: debugR,
          loading: false,
          plotPath: output.plotPath,
          txtPath: output.txtPath,
        });
        setRPlot(output.plotPath);
      }
    } catch (err) {
      dispatchError(err);
      dispatchExpPrevalence({ loading: false });
    }
  }

  async function setRPlot(plotPath) {
    if (plotPath) {
      try {
        const response = await fetch(`${rootURL}results/${plotPath}`);
        if (!response.ok) {
          // console.log(await response.json());
        } else {
          const pic = await response.blob();
          const objectURL = URL.createObjectURL(pic);

          if (plotURL) URL.revokeObjectURL(plotURL);
          dispatchExpPrevalence({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchExpPrevalence({ err: true, plotURL: '' });
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

    dispatchExpPrevalence({
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
                id="prevalenceStudy"
                label="Study"
                value={study}
                options={studyOptions}
                onChange={handleStudy}
              />
            </Col>
            <Col sm="2">
              <Select
                id="prevalenceStrategy"
                label="Experimental Strategy"
                value={strategy}
                options={strategyOptions}
                onChange={(strategy) =>
                  dispatchExpPrevalence({ strategy: strategy })
                }
              />
            </Col>
            <Col sm="3">
              <Select
                id="prevalenceRefSet"
                label="Reference Signature Set"
                value={refSignatureSet}
                options={refSignatureSetOptions}
                onChange={(set) =>
                  dispatchExpPrevalence({ refSignatureSet: set })
                }
              />
            </Col>
            <Col sm="2">
              <Group controlId="prevalenceGenomeSize">
                <Label>Genome Size</Label>
                <Control
                  value={genomeSize}
                  onChange={(e) => {
                    dispatchExpPrevalence({
                      genomeSize: e.target.value,
                    });
                  }}
                ></Control>
                {/* <Text className="text-muted">(Ex. NCG>NTG)</Text> */}
              </Group>
            </Col>
            <Col sm="2">
              <Group controlId="prevalenceMutations">
                <Label>Minimal Number Mutations within in Each Signature</Label>
                <Control
                  value={mutation}
                  onChange={(e) => {
                    dispatchExpPrevalence({
                      mutation: e.target.value,
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
                  calculateR('cosineSimilarity', {
                    study: study,
                    strategy: strategy,
                    refSignatureSet: refSignatureSet,
                    genomeSize: parseFloat(genomeSize),
                    mutation: parseFloat(mutation),
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
