import React, { useState } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Plot from '../../controls/plot/plot';
import CustomSelect from '../../controls/select/select';
import Description from '../../controls/description/description';
import KataegisTable from './kataegisTable';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../services/store/visualization';
import { actions as modalActions } from '../../../services/store/modal';
import { NavHashLink } from 'react-router-hash-link';

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
            kataegisData: data ? JSON.parse(data) : [],
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
        <div className="bg-white border rounded">
          <Description
            className="p-3"
            less="This analysis identifies the kataegis events from a VCF file input."
            more={
              <>
                <span>
                  Kataegis is a localized substitution hypermutation event,
                  often characterized by clusters of C>T and/or C>G mutations,
                  commonly at TpCpN trinucleotides (APOBEC mutations). Click{' '}
                  <NavHashLink to="/faq#kataegis">here</NavHashLink> for
                  additional information about kataegis.
                </span>
                <p className="mt-3">
                  To identify kataegis, input the [Minimum Number of Mutations]
                  required for kataegis, the [Maximum Distance] between one
                  mutation and the next within a given cluster of mutations
                  being considered for kataegis, and a [Chromosome] to be
                  highlighted in the rainfall plot. By default, all chromosomes
                  will be shown for the kataegis identification.
                </p>
              </>
            }
          />

          <hr />
          <Form noValidate className="p-3">
            <LoadingOverlay active={loading} />
            <Row>
              <Col lg="auto">
                <CustomSelect
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
                <CustomSelect
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
              <>
                <hr /> <p className="p-3 text-danger">{err}</p>
              </>
            )}
            {plotPath && (
              <>
                <hr />
                <Plot
                  className="p-3"
                  downloadName={plotPath.split('/').slice(-1)[0]}
                  plotPath={'api/results/' + plotPath}
                  txtPath={txtPath ? `api/results/${txtPath}` : null}
                />
                <p className="p-3">
                  The rainfall plot illustrates the kataegis events identified
                  given the input parameters. The x-axis is the position of the
                  mutation in the genome separated by chromosome. The y-axis is
                  the inter-mutation distance (bp, log10) from one mutation to
                  the next. The colors represent different mutation types. The
                  arrow will highlight the kataegis region when you select the
                  “Highlight” option.
                </p>
              </>
            )}
            {kataegisData.length > 0 && (
              <>
                <p className="p-3">
                  The table below is a summary of the kataegis identification
                  based on the input parameters. The table can be filtered based
                  on any of the columns by entering an appropriate value.
                </p>
                <KataegisTable />
              </>
            )}
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
