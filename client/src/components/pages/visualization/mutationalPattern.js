import React, { useEffect, useState } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Plot from '../../controls/plot/plot';
import Debug from '../../controls/debug/debug';
import { value2d, filter2d } from '../../../services/utils';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../services/store/visualization';
import { actions as modalActions } from '../../../services/store/modal';

const actions = { ...visualizationActions, ...modalActions };
const { Group, Label, Control } = Form;

export default function MutationalPattern({ submitR }) {
  const dispatch = useDispatch();
  const visualization = useSelector((state) => state.visualization);
  const mergeMPEA = (state) =>
    dispatch(actions.mergeVisualization({ mutationalPattern: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const {
    source,
    study,
    cancerType,
    pubExperimentalStrategy,
  } = visualization.visualize;
  const { projectID, matrixList } = visualization.results;
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
  } = visualization.mutationalPattern;

  const [invalidProportion, setProportion] = useState(false);
  const [invalidPattern, setPattern] = useState(false);

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
        mergeMPEA({ barPath: barchart });
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
        mergeMPEA({
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
    mergeMPEA({ loading: true });

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
            mergeMPEA({
              plotURL: objectURL,
            });
          } else {
            if (barURL) URL.revokeObjectURL(barURL);
            mergeMPEA({
              barURL: objectURL,
            });
          }
        }
      } catch (err) {
        mergeError(err.message);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      mergeMPEA({ err: true, plotURL: '' });
    }
    mergeMPEA({ loading: false });
  }

  function clearPlot(type) {
    if (type == 'context') {
      URL.revokeObjectURL(plotURL);
      mergeMPEA({ plotURL: '' });
    } else if (type == 'barchart') {
      URL.revokeObjectURL(barURL);
      mergeMPEA({ barURL: '' });
    }
  }

  async function calculateR(fn, args) {
    mergeMPEA({
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
        mergeMPEA({ debugR: err, loading: false });
      } else {
        const { debugR, output } = await response.json();
        if (Object.keys(output).length) {
          mergeMPEA({
            debugR: debugR,
            plotPath: output.plotPath,
            barPath: output.barPath,
            txtPath: output.txtPath,
          });
        } else {
          mergeMPEA({
            debugR: debugR,
            err: true,
            loading: false,
          });
        }
      }
    } catch (err) {
      mergeError(err.message);
      mergeMPEA({ loading: false });
    }
  }

  const plots = (
    <div>
      {err && (
        <div>
          <hr />
          <p className="p-3">An error has occured. Please verify your input.</p>
        </div>
      )}

      <div id="barchart">
        {barPath && (
          <>
            <hr />
            <Plot
              className="p-3"
              downloadName={barPath.split('/').slice(-1)[0]}
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
              downloadName={plotPath.split('/').slice(-1)[0]}
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
      <div className="bg-white border rounded">
        <Form noValidate className="p-3">
          <LoadingOverlay active={loading} />
          <Row>
            <Col lg="2">
              <Group
                controlId="minimum"
                title="Minimal Proportion mutations within Each Mutational Pattern"
              >
                <Label>Minimal Proportion</Label>
                <Control
                  value={proportion}
                  placeholder="Ex. 0.8"
                  onChange={(e) => {
                    mergeMPEA({
                      proportion: e.target.value,
                    });
                  }}
                  isInvalid={invalidProportion}
                />
                <Form.Control.Feedback type="invalid">
                  Enter a valid proportion between 0 and 1
                </Form.Control.Feedback>
              </Group>
            </Col>
            <Col lg="2">
              <Group controlId="pattern">
                <Label>Mutational Pattern</Label>
                <Control
                  value={pattern}
                  placeholder="Ex. NCG>NTG"
                  onChange={(e) => {
                    mergeMPEA({
                      pattern: e.target.value,
                    });
                  }}
                  isInvalid={invalidPattern}
                />
                <Form.Control.Feedback type="invalid">
                  Enter a valid pattern
                </Form.Control.Feedback>
              </Group>
            </Col>
            <Col />
            <Col lg="2" className="d-flex justify-content-end">
              <Button
                className="mt-auto mb-3"
                variant="primary"
                onClick={() => {
                  if (!pattern) setPattern(true);
                  else setPattern(false);
                  if (isNaN(proportion) || !proportion) setProportion(true);
                  else setProportion(false);

                  if (pattern && proportion && !isNaN(proportion)) {
                    if (source == 'user') {
                      calculateR('mutationalPattern', {
                        matrixFile: value2d(
                          filter2d(['SBS', '96'], matrixList.data)[0],
                          'Path',
                          matrixList.columns
                        ),
                        proportion: parseFloat(proportion),
                        pattern: pattern,
                      });
                    } else if (source == 'public') {
                      calculateR('mutationalPatternPublic', {
                        study: study,
                        cancerType: cancerType,
                        experimentalStrategy: pubExperimentalStrategy,
                        proportion: parseFloat(proportion),
                        pattern: pattern,
                      });
                    }
                  }
                }}
              >
                Calculate
              </Button>
            </Col>
          </Row>
        </Form>
        {plots}
      </div>

      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
