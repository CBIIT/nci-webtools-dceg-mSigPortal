import React, { useEffect, useState } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Plot from '../../controls/plot/plot';
import Debug from '../../controls/debug/debug';
import Select from '../../controls/select/select';
import KataegisTable from './kataegisTable';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../services/store/visualization';
import { actions as modalActions } from '../../../services/store/modal';

const actions = { ...visualizationActions, ...modalActions };
const { Group, Check, Label, Control } = Form;

export default function Kataegis({ submitR }) {
  const dispatch = useDispatch();
  const visualization = useSelector((state) => state.visualization);
  const mergeKataegis = (state) =>
    dispatch(actions.mergeVisualization({ kataegis: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const { source, inputFormat } = visualization.visualize;
  const { projectID } = visualization.results;
  const {
    sample,
    sampleOptions,
    highlight,
    min,
    max,
    chromosome,
    txtPath,
    plotPath,
    plotURL,
    err,
    kataegisData,
    debugR,
    loading,
  } = visualization.kataegis;

  const [invalidMin, setMin] = useState(false);
  const [invalidMax, setMax] = useState(false);

  useEffect(() => {
    if (plotPath) setRPlot(plotPath);
    else clearPlot();
  }, [plotPath]);

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
          mergeKataegis({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        mergeError(err.message);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      mergeKataegis({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      mergeKataegis({ plotPath: '', plotURL: '' });
    }
  }

  async function calculateR() {
    mergeKataegis({
      loading: true,
      err: false,
      debugR: '',
      plotPath: '',
      kataegisData: [],
    });

    try {
      const response = await submitR('kataegis', {
        sample: sample,
        highlight: highlight,
        min: parseInt(min),
        max: parseInt(max),
        chromosome: chromosome,
      });
      if (!response.ok) {
        const err = await response.json();
        mergeKataegis({ debugR: err });

        mergeKataegis({ loading: false });
      } else {
        const { debugR, output, errors, data } = await response.json();
        if (Object.keys(output).length) {
          mergeKataegis({
            debugR: debugR,
            plotPath: output.plotPath,
            txtPath: output.txtPath,
            kataegisData: JSON.parse(data),
            loading: false,
          });
        } else {
          mergeKataegis({
            debugR: debugR,
            plotPath: '',
            txtPath: '',
            err: errors,
            kataegisData: [],
            loading: false,
          });
        }
      }
    } catch (err) {
      mergeError(err.message);
      mergeKataegis({ loading: false });
    }
  }

  return (
    <div>
      {source == 'user' && inputFormat == 'vcf' ? (
        <div className="bg-white border rounded">
          <Form noValidate className="p-3">
            <LoadingOverlay active={loading} />
            <Row>
              <Col lg="3">
                <Select
                  id="kataegisSamples"
                  label="Sample Name"
                  value={sample}
                  options={sampleOptions}
                  onChange={(sample) =>
                    mergeKataegis({
                      sample: sample,
                    })
                  }
                />
              </Col>
              <Col lg="1">
                <Group controlId="toggleHighlight">
                  <Check
                    type="checkbox"
                    label="Highlight"
                    value={highlight}
                    checked={highlight}
                    onChange={() => {
                      mergeKataegis({ highlight: !highlight });
                    }}
                  />
                </Group>
              </Col>
              <Col lg="2">
                <Group controlId="kataegisMin">
                  <Label>Minimum Number of Mutations</Label>
                  <Control
                    value={min}
                    placeholder=""
                    onChange={(e) => {
                      mergeKataegis({
                        min: e.target.value,
                      });
                    }}
                    isInvalid={invalidMin}
                  />
                  <Form.Control.Feedback type="invalid">
                    Enter a numeric minimum value
                  </Form.Control.Feedback>
                </Group>
              </Col>
              <Col lg="2">
                <Group controlId="kataegisMax">
                  <Label>Maximum Distance</Label>
                  <Control
                    value={max}
                    placeholder=""
                    onChange={(e) => {
                      mergeKataegis({
                        max: e.target.value,
                      });
                    }}
                    isInvalid={invalidMax}
                  />
                  <Form.Control.Feedback type="invalid">
                    Enter a numeric maximum value
                  </Form.Control.Feedback>
                </Group>
              </Col>
              <Col lg="2">
                <Select
                  id="kataegisChromosome"
                  label="Chromosome"
                  value={chromosome}
                  options={[
                    'None',
                    ...[...Array(23).keys()].slice(1).map((chr) => 'chr' + chr),
                    'chrX',
                    'chrY',
                  ]}
                  onChange={(chromosome) =>
                    mergeKataegis({
                      chromosome: chromosome,
                    })
                  }
                />
              </Col>
              <Col lg="2" className="d-flex justify-content-end">
                <Button
                  className="mt-auto mb-3"
                  variant="primary"
                  onClick={() => {
                    if (!min || isNaN(min)) setMin(true);
                    else setMin(false);
                    if (!max || isNaN(max)) setMax(true);
                    else setMax(false);

                    if (min && max && !isNaN(min) && !isNaN(max)) calculateR();
                  }}
                >
                  Calculate
                </Button>
              </Col>
            </Row>
          </Form>

          <div id="kataegisPlot">
            {err && (
              <div>
                <hr />
                <p className="p-3 text-danger">{err}</p>
              </div>
            )}
            <div style={{ display: plotURL ? 'block' : 'none' }}>
              <hr />
              <Plot
                className="p-3"
                downloadName={plotPath.split('/').slice(-1)[0]}
                plotURL={plotURL}
                txtPath={txtPath ? projectID + txtPath : null}
              />
              {kataegisData.length > 0 && <KataegisTable />}
            </div>
          </div>
        </div>
      ) : (
        <div className="bg-white border rounded p-4">
          <p>
            Kataegis Identification is only available for Data Source: User and
            requires VCF input files
          </p>
        </div>
      )}
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
