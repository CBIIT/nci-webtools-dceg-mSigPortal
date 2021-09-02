import React, { useState } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Plot from '../../controls/plot/plot';
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
    projectID,
    matrixList,
  } = visualization.state;

  const {
    proportion,
    pattern,
    txtPath,
    plotPath,
    barPath,
    err,
    debugR,
    loading,
  } = visualization.mutationalPattern;

  const [tmpProportion, setProportion] = useState(proportion);
  const [invalidProportion, setProportionInvalid] = useState(false);
  const [invalidPattern, setPatternInvalid] = useState(false);

  async function calculateR(fn, args) {
    mergeMPEA({
      loading: true,
      err: false,
      debugR: '',
      plotPath: '',
      txtPath: '',
      barPath: '',
      proportion: tmpProportion,
    });

    try {
      const response = await submitR(fn, args);
      if (!response.ok) {
        const err = await response.json();
        mergeMPEA({ debugR: err });
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
          });
        }
      }
    } catch (err) {
      mergeError(err.message);
    }
    mergeMPEA({ loading: false });
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
              plotPath={'api/results/' + projectID + barPath}
            />
            <p className="p-3">
              This plot illustrates the frequency by count of each mutational
              pattern in the given study and cancer type or input dataset. The
              y-axis is the frequency of each mutational pattern across all
              samples, and the x-axis includes each of the mutational patterns
              present in the study and cancer type that meet the criteria for
              the minimal proportion of mutations within each mutational
              pattern.
            </p>
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
              plotPath={'api/results/' + projectID + plotPath}
              txtPath={projectID + txtPath}
              title="Proportion of Mutational Pattern Context Compared to Other Contexts with the same SBS Mutation"
            />
            <p className="p-3">
              This plot illustrates the mutational pattern context entered
              compared to other contexts with the same SBS mutation for each
              sample. On the y-axis is the other contexts, and on the x-axis is
              the specific mutational pattern context input for the enrichment
              analysis. For some studies including multiple cancer types (such
              as TCGA PanCancer), different colors will be used.
            </p>
          </>
        )}
      </div>
    </div>
  );

  return (
    <div>
      <div className="bg-white border rounded">
        <div className="p-3">
          <p>
            This page allows you to conduct a mutational pattern enrichment
            analysis. This type of analysis aims to determine frequency and
            enrichment of different types of mutational patterns. For more
            information about mutational pattern enrichment, click{' '}
            <a href="#faq">here</a>. Below are explanations of the inputs needed
            for the analysis.
          </p>
          <p>
            Minimal Proportion: For “Frequency of Mutational Pattern” plot, set
            the minimal proportion of mutational patterns identified in all
            samples from selected or input study. A slightly high proportion
            like 0.5 are suggested.{' '}
          </p>
          <p>
            Mutational Pattern: For the second enrichment plot, select the
            mutational patten to identify the enrichment of specific mutation
            context as suggested from “Frequency of Mutational Patten”. The
            mutational pattern supports common nucleotide symbols.{' '}
          </p>
        </div>
        <hr />
        <Form noValidate className="p-3">
          <LoadingOverlay active={loading} />
          <Row>
            <Col lg="auto">
              <Group
                controlId="minimum"
                title="Minimal Proportion mutations within Each Mutational Pattern"
              >
                <Label>Minimal Proportion (0-1)</Label>
                <Control
                  value={tmpProportion}
                  placeholder="Ex. 0.8"
                  onChange={(e) => setProportion(e.target.value)}
                  isInvalid={invalidProportion}
                />
                <Form.Control.Feedback type="invalid">
                  Please input a value form 0 to 1 for Minimal Proportion
                </Form.Control.Feedback>
              </Group>
            </Col>
            <Col lg="auto">
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
            <Col lg="2" className="d-flex">
              <Button
                className="ml-auto mb-auto"
                variant="primary"
                onClick={() => {
                  if (!pattern) setPatternInvalid(true);
                  else setPatternInvalid(false);
                  if (
                    isNaN(tmpProportion) ||
                    !tmpProportion ||
                    tmpProportion < 0 ||
                    tmpProportion > 1
                  )
                    setProportionInvalid(true);
                  else setProportionInvalid(false);

                  if (
                    pattern &&
                    tmpProportion &&
                    !isNaN(tmpProportion) &&
                    tmpProportion >= 0 &&
                    tmpProportion <= 1
                  ) {
                    if (source == 'user') {
                      calculateR('mutationalPattern', {
                        matrixFile: matrixList.filter(
                          (row) =>
                            row.Profile_Type == 'SBS' && row.Matrix_Size == '96'
                        )[0].Path,
                        proportion: parseFloat(tmpProportion),
                        pattern: pattern,
                      });
                    } else if (source == 'public') {
                      calculateR('mutationalPatternPublic', {
                        study: study,
                        cancerType: cancerType,
                        experimentalStrategy: pubExperimentalStrategy,
                        proportion: parseFloat(tmpProportion),
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
