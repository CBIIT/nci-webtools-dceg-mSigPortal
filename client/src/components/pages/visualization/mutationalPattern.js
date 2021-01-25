import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchMutationalPattern,
} from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Plot from '../../controls/plot/plot';
import Debug from '../../controls/debug/debug';

const { Group, Label, Control } = Form;

export default function MutationalPattern({ submitR }) {
  const { source, study, cancerType, pubExperimentalStrategy } = useSelector(
    (state) => state.visualize
  );
  const { projectID, matrixList } = useSelector(
    (state) => state.visualizeResults
  );
  const {
    proportion,
    pattern,
    txtPath,
    plotPath,
    plotURL,
    barPath,
    barURL,
    err,
    debugR,
    loading,
  } = useSelector((state) => state.mutationalPattern);

  // load plots if they exist - used for precalculated examples
  useEffect(() => {
    const checkPlot = async () => {
      const barchart =
        source == 'user'
          ? `/results/mutationalPattern/barchart.svg`
          : `/results/mutationalPatternPublic/barchart.svg`;
      let check = await fetch(`api/results/${projectID}${barchart}`, {
        method: 'HEAD',
        cache: 'no-cache',
      });
      if (check.status === 200) {
        dispatchMutationalPattern({ barPath: barchart });
      }

      const mpea =
        source == 'user'
          ? `/results/mutationalPattern/mpea.svg`
          : `/results/mutationalPatternPublic/mpea.svg`;
      check = await fetch(`api/results/${projectID}${mpea}`, {
        method: 'HEAD',
        cache: 'no-cache',
      });
      if (check.status === 200) {
        dispatchMutationalPattern({
          plotPath: mpea,
          txtPath: `/results/mutationalPatternPublic/mpea.txt`,
        });
      }
    };

    if (projectID) checkPlot();
  }, [projectID]);

  useEffect(() => {
    plotPath ? setRPlot(plotPath, 'context') : clearPlot('context');
    barPath ? setRPlot(barPath, 'barchart') : clearPlot('barchart');
  }, [plotPath, barPath]);

  async function setRPlot(plotPath, type) {
    dispatchMutationalPattern({ loading: true });

    if (plotPath) {
      try {
        const response = await fetch(`api/results/${projectID}${plotPath}`);
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
    dispatchMutationalPattern({ loading: false });
  }

  function clearPlot(type) {
    if (type == 'context') {
      URL.revokeObjectURL(plotURL);
      dispatchMutationalPattern({ plotURL: '' });
    } else if (type == 'barchart') {
      URL.revokeObjectURL(barURL);
      dispatchMutationalPattern({ barURL: '' });
    }
  }

  async function calculateR(fn, args) {
    dispatchMutationalPattern({
      loading: true,
      err: false,
      debugR: '',
      plotPath: '',
      txtPath: '',
      barPath: '',
    });

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();
        dispatchMutationalPattern({ debugR: err, loading: false });
      } else {
        const { debugR, output } = await response.json();
        if (Object.keys(output).length) {
          dispatchMutationalPattern({
            debugR: debugR,
            plotPath: output.plotPath,
            barPath: output.barPath,
            txtPath: output.txtPath,
          });
        } else {
          dispatchMutationalPattern({
            debugR: debugR,
            err: true,
            loading: false,
          });
        }
      }
    } catch (err) {
      dispatchError(err);
      dispatchMutationalPattern({ loading: false });
    }
  }

  const plots = (
    <div>
      {err && (
        <div className="p-3">
          <p>An error has occured. Please verify your input.</p>
        </div>
      )}

      <div id="barchart">
        {barPath && (
          <>
            <hr />
            <Plot
              className="p-3"
              plotName={barPath.split('/').slice(-1)[0]}
              plotURL={barURL}
            />
          </>
        )}
        {plotPath && !barPath && (
          <div className="p-3">
            <p>Frequency of Mutational Pattern</p>
            <p>
              No mutational pattern with proportion of mutations large than{' '}
              {proportion}
            </p>
          </div>
        )}
      </div>
      <div id="context">
        {plotPath && (
          <>
            <hr />
            <Plot
              className="p-3"
              plotName={plotPath.split('/').slice(-1)[0]}
              plotURL={plotURL}
              txtPath={projectID + txtPath}
            />
          </>
        )}
      </div>
    </div>
  );

  return (
    <div>
      {source == 'user' ? (
        <div className="bg-white border rounded">
          <Form noValidate className="p-3">
            <LoadingOverlay active={loading} />
            <Row>
              <Col lg="4">
                <Group controlId="minimum">
                  <Label>
                    Minimal Proportion mutations within Each Mutational Pattern
                  </Label>
                  <Control
                    value={proportion}
                    placeholder="Ex. 0.8"
                    onChange={(e) => {
                      dispatchMutationalPattern({
                        proportion: e.target.value,
                      });
                    }}
                    isInvalid={!proportion || isNaN(proportion)}
                  />
                  <Form.Control.Feedback type="invalid">
                    Enter a valid proportion between 0 and 1
                  </Form.Control.Feedback>
                </Group>
              </Col>
              <Col lg="3">
                <Group controlId="pattern">
                  <Label>Mutational Pattern</Label>
                  <Control
                    value={pattern}
                    placeholder="Ex. NCG>NTG"
                    onChange={(e) => {
                      dispatchMutationalPattern({
                        pattern: e.target.value,
                      });
                    }}
                    isInvalid={!pattern}
                  />
                  <Form.Control.Feedback type="invalid">
                    Enter a valid pattern
                  </Form.Control.Feedback>
                </Group>
              </Col>
              <Col lg="3" />
              <Col lg="2" className="d-flex justify-content-end">
                <Button
                  className="mt-auto mb-3"
                  variant="primary"
                  onClick={() => {
                    calculateR('mutationalPattern', {
                      matrixFile: matrixList.filter(
                        (path) =>
                          path.Profile_Type == 'SBS' && path.Matrix_Size == '96'
                      )[0].Path,
                      proportion: parseFloat(proportion),
                      pattern: pattern,
                    });
                  }}
                  disabled={!pattern || isNaN(proportion) || !proportion}
                >
                  Calculate
                </Button>
              </Col>
            </Row>
          </Form>
          {plots}
        </div>
      ) : (
        <div className="bg-white border rounded">
          <Form className="p-3">
            <LoadingOverlay active={loading} />
            <Row className="justify-content-center">
              <Col lg="4">
                <Label>
                  Minimal Proportion mutations within Each Mutational Pattern
                </Label>
                <Control
                  value={proportion}
                  onChange={(e) => {
                    dispatchMutationalPattern({
                      proportion: e.target.value,
                    });
                  }}
                  isInvalid={!proportion || isNaN(proportion)}
                />
                <Form.Control.Feedback type="invalid">
                  Enter a valid proportion between 0 and 1
                </Form.Control.Feedback>
              </Col>
              <Col lg="3">
                <Label>Mutational Pattern</Label>
                <Control
                  value={pattern}
                  onChange={(e) => {
                    dispatchMutationalPattern({
                      pattern: e.target.value,
                    });
                  }}
                  isInvalid={!pattern}
                />
                <Form.Control.Feedback type="invalid">
                  Enter a valid pattern
                </Form.Control.Feedback>
              </Col>
              <Col lg="3" />
              <Col lg="2" className="d-flex justify-content-end">
                <Button
                  className="mt-auto mb-3"
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
                  disabled={!pattern || isNaN(proportion) || !proportion}
                >
                  Calculate
                </Button>
              </Col>
            </Row>
          </Form>
          {plots}
        </div>
      )}
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
