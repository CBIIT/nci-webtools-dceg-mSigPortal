import React, { useState } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Plot from '../../controls/plot/plot';
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

  const { source, inputFormat, projectID } = visualization.state;

  const {
    sample,
    sampleOptions,
    highlight,
    min,
    max,
    chromosome,
    txtPath,
    plotPath,
    err,
    kataegisData,
    loading,
  } = visualization.kataegis;

  const [invalidMin, setMin] = useState(false);
  const [invalidMax, setMax] = useState(false);

  async function calculateR() {
    mergeKataegis({
      loading: true,
      err: false,
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
        mergeKataegis({ err: err });
        mergeKataegis({ loading: false });
      } else {
        const { output } = await response.json();
        const { plotPath, txtPath, data, error, uncaughtError } = output;
        if (plotPath) {
          mergeKataegis({
            plotPath: plotPath,
            txtPath: txtPath,
            kataegisData: JSON.parse(data),
            loading: false,
          });
        } else {
          mergeKataegis({
            plotPath: '',
            txtPath: '',
            err: error || uncaughtError,
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
        <div className="bg-white border rounded" style={{ minHeight: '500px' }}>
          <div className="p-3">
            <p>
              Below you can investigate instances of kataegis using a VCF file.
              Kataegis is localized substitution hypermutation, often
              characterized by clusters of C>T and/or C>G mutations, commonly at
              TpCpN trinucleotides. Click <a href="#faq">here</a> for additional
              information about kataegis. Simply input the “Minimum Number of
              Mutations” to qualify was kataegis, the “Maximum Distance” between
              one mutation and the next within a given cluster of mutations
              being considered for kataegis, and a “Chromosome” to be included
              in the Kataegis results figure. If you would like for all
              chromosomes to be output in the results table, leave the
              “Chromosome” input as ‘None’.
            </p>
          </div>
          <hr />
          <Form noValidate className="p-3">
            <LoadingOverlay active={loading} />
            <Row>
              <Col lg="auto">
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
              <Col lg="auto">
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
              <Col lg="auto">
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
              <Col lg="auto">
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
              <Col lg="auto">
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
              <Col lg="auto" className="d-flex">
                <Button
                  className="ml-auto mb-auto"
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
          <hr />

          <div id="kataegisPlot">
            {err && <p className="p-3 text-danger">{err}</p>}
            <div style={{ display: plotPath ? 'block' : 'none' }}>
              <Plot
                className="p-3"
                downloadName={plotPath.split('/').slice(-1)[0]}
                plotPath={'api/results/' + plotPath}
                txtPath={txtPath ? `api/results/${txtPath}` : null}
              />
              <p className="p-3">
                The rainfall plot illustrates the kataegis identified given the
                input parameters. Along the y-axis is the inter-mutation
                distance (bp, log10) from one mutation to the next. On the
                x-axis is the position of the mutation in the genome by
                chromosome. The colors represent the different mutation types.
                The arrow will highlight the kataegis region when you select the
                “Highlight”.
              </p>
              {kataegisData.length > 0 && (
                <>
                  <p className="p-3">
                    The table below is a summary of the kataegis identification
                    for the input file. The table can be filtered based on any
                    of the columns by entering an appropriate value for the
                    given column in the box above each column.{' '}
                  </p>
                  <KataegisTable />
                </>
              )}
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
    </div>
  );
}
