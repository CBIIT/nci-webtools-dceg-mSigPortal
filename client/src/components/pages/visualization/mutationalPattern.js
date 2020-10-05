import React from 'react';
import { Form, Row, Col, Button, Accordion, Card } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faPlus, faMinus } from '@fortawesome/free-solid-svg-icons';
import {
  dispatchError,
  dispatchMutationalPattern,
} from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Plot from '../../controls/plot/plot';
import Debug from '../../controls/debug/debug';

const { Group, Label, Control, Text } = Form;
const { Header, Body } = Card;
const { Toggle, Collapse } = Accordion;

export default function MutationalPattern({ submitR }) {
  const { source, study, cancerType, pubExperimentalStrategy } = useSelector(
    (state) => state.visualize
  );
  const { matrixList } = useSelector((state) => state.visualizeResults);
  const rootURL = window.location.pathname;
  const {
    proportion,
    pattern,
    txtPath,
    plotPath,
    plotURL,
    barPath,
    barURL,
    display,
    err,
    debugR,
    displayDebug,
    submitOverlay,
  } = useSelector((state) => state.mutationalPattern);

  async function setRPlot(plotPath, type) {
    dispatchMutationalPattern({ submitOverlay: true });

    if (plotPath) {
      try {
        const response = await fetch(`${rootURL}results/${plotPath}`);
        if (!response.ok) {
          // console.log(await response.json());
        } else {
          const pic = await response.blob();
          const objectURL = URL.createObjectURL(pic);

          if (type == 'context') {
            if (plotURL) URL.revokeObjectURL(plotURL);
            dispatchMutationalPattern({
              plotURL: objectURL,
            });
          } else {
            if (barURL) URL.revokeObjectURL(barURL);
            dispatchMutationalPattern({
              barURL: objectURL,
            });
          }
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchMutationalPattern({ err: true, plotURL: '' });
    }
    dispatchMutationalPattern({ submitOverlay: false });
  }

  async function calculateR(fn, args) {
    dispatchMutationalPattern({
      submitOverlay: true,
      err: false,
      debugR: '',
    });

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();
        dispatchMutationalPattern({ debugR: err });

        dispatchMutationalPattern({ submitOverlay: false });
      } else {
        const { debugR, output } = await response.json();
        if (output.plotPath) {
          dispatchMutationalPattern({
            debugR: debugR,
            plotPath: output.plotPath,
            barPath: output.barPath,
            txtPath: output.txtPath,
          });
          setRPlot(output.plotPath, 'context');
          if (output.barPath) setRPlot(output.barPath, 'barchart');
        } else {
          dispatchMutationalPattern({
            debugR: debugR,
            err: true,
            submitOverlay: false,
          });
        }
      }
    } catch (err) {
      dispatchError(err);
      dispatchMutationalPattern({ submitOverlay: false });
    }
  }

  const plots = (
    <>
      <div id="barchart">
        <div style={{ display: err ? 'block' : 'none' }}>
          <h4>Bar Chart</h4>
          <p>An error has occured. Check the debug section for more info.</p>
        </div>
        {plotURL.length > 0 &&
          (barURL.length > 0 ? (
            <Plot plotName={barPath.split('/').slice(-1)[0]} plotURL={barURL} />
          ) : (
            <div>
              <h4>Proportion</h4>
              <p>
                No mutational pattern with proportion of mutations large than{' '}
                {proportion}
              </p>
            </div>
          ))}
      </div>
      <div id="context">
        <div style={{ display: err ? 'block' : 'none' }}>
          <h4>Context Plot</h4>
          <p>An error has occured. Check the debug section for more info.</p>
        </div>
        <div style={{ display: plotURL ? 'block' : 'none' }}>
          <Plot
            plotName={plotPath.split('/').slice(-1)[0]}
            plotURL={plotURL}
            txtPath={txtPath}
          />
        </div>
      </div>
    </>
  );

  return (
    <div>
      {source == 'user' && (
        <Accordion defaultActiveKey="0">
          <Card>
            <Toggle
              className="font-weight-bold"
              as={Header}
              eventKey="0"
              onClick={() =>
                dispatchMutationalPattern({
                  display: !display,
                })
              }
            >
              {display == true ? (
                <FontAwesomeIcon icon={faMinus} />
              ) : (
                <FontAwesomeIcon icon={faPlus} />
              )}{' '}
              Motif Analysis
            </Toggle>
            <Collapse eventKey="0">
              <Body>
                <Form>
                  <LoadingOverlay active={submitOverlay} />
                  <div>
                    <Row className="justify-content-center">
                      <Col sm="5">
                        <Group controlId="minimum">
                          <Label>
                            Minimal Proportion mutations within Each Mutational
                            Pattern
                          </Label>
                          <Control
                            value={proportion}
                            onChange={(e) => {
                              dispatchMutationalPattern({
                                proportion: e.target.value,
                              });
                            }}
                          ></Control>{' '}
                          <Text className="text-muted">(Ex. 0.8)</Text>
                        </Group>
                      </Col>
                      <Col sm="5">
                        <Group controlId="pattern">
                          <Label>Mutational Pattern</Label>
                          <Control
                            value={pattern}
                            onChange={(e) => {
                              dispatchMutationalPattern({
                                pattern: e.target.value,
                              });
                            }}
                          ></Control>
                          <Text className="text-muted">(Ex. NCG>NTG)</Text>
                        </Group>
                      </Col>

                      <Col sm="1" className="m-auto">
                        <Button
                          variant="primary"
                          onClick={() => {
                            calculateR('mutationalPattern', {
                              matrixFile: matrixList.filter(
                                (path) =>
                                  path.Profile_Type == 'SBS' &&
                                  path.Matrix_Size == '96'
                              )[0].Path,
                              proportion: parseFloat(proportion),
                              pattern: pattern,
                            });
                          }}
                        >
                          Calculate
                        </Button>
                      </Col>
                    </Row>{' '}
                    {plots}
                  </div>
                </Form>
              </Body>
            </Collapse>
          </Card>
        </Accordion>
      )}

      {source == 'public' && (
        <Accordion defaultActiveKey="1">
          <Card>
            <Toggle
              className="font-weight-bold"
              as={Header}
              eventKey="1"
              onClick={() =>
                dispatchMutationalPattern({
                  display: !display,
                })
              }
            >
              {display == true ? (
                <FontAwesomeIcon icon={faMinus} />
              ) : (
                <FontAwesomeIcon icon={faPlus} />
              )}{' '}
              Motif Analysis
            </Toggle>
            <Collapse eventKey="1">
              <Body>
                <Form>
                  <LoadingOverlay active={submitOverlay} />
                  <div>
                    <Row className="justify-content-center">
                      <Col sm="5">
                        <Label>
                          Minimal Proportion mutations within Each Mutational
                          Pattern
                        </Label>
                        <Control
                          value={proportion}
                          onChange={(e) => {
                            dispatchMutationalPattern({
                              proportion: e.target.value,
                            });
                          }}
                        ></Control>{' '}
                        <Text className="text-muted">(Ex. 0.8)</Text>
                      </Col>
                      <Col sm="5">
                        <Label>Mutational Pattern</Label>
                        <Control
                          value={pattern}
                          onChange={(e) => {
                            dispatchMutationalPattern({
                              pattern: e.target.value,
                            });
                          }}
                        ></Control>
                        <Text className="text-muted">(Ex. NCG>NTG)</Text>
                      </Col>

                      <Col sm="1" className="m-auto">
                        <Button
                          variant="primary"
                          onClick={() => {
                            calculateR('mutationalPatternPublic', {
                              study: study,
                              cancerType: cancerType,
                              experimentalStrategy: pubExperimentalStrategy,
                              proportion: parseFloat(proportion),
                              pattern: pattern,
                            });
                          }}
                        >
                          Calculate
                        </Button>
                      </Col>
                    </Row>
                    {plots}
                  </div>
                </Form>
              </Body>
            </Collapse>
          </Card>
        </Accordion>
      )}
      <Debug msg={debugR} />
    </div>
  );
}
