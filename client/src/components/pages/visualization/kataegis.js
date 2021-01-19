import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { dispatchError, dispatchKataegis } from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Plot from '../../controls/plot/plot';
import Debug from '../../controls/debug/debug';
import Select from '../../controls/select/select';
import KataegisTable from './kataegisTable';

const { Group, Check, Label, Control } = Form;

export default function Kataegis({ submitR }) {
  const { source, inputFormat } = useSelector((state) => state.visualize);
  const { projectID } = useSelector((state) => state.visualizeResults);
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
    debugR,
    loading,
  } = useSelector((state) => state.kataegis);

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
          dispatchKataegis({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchKataegis({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      dispatchKataegis({ plotPath: '', plotURL: '' });
    }
  }

  async function calculateR() {
    dispatchKataegis({
      loading: true,
      err: false,
      debugR: '',
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
        dispatchKataegis({ debugR: err });

        dispatchKataegis({ loading: false });
      } else {
        const { debugR, output, errors, data } = await response.json();
        if (Object.keys(output).length) {
          dispatchKataegis({
            debugR: debugR,
            plotPath: output.plotPath,
            txtPath: output.txtPath,
            kataegisData: JSON.parse(data),
            loading: false,
          });
        } else {
          dispatchKataegis({
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
      dispatchError(err);
      dispatchKataegis({ loading: false });
    }
  }

  return (
    <div>
      {source == 'user' && inputFormat == 'vcf' ? (
        <div className="bg-white border rounded">
          <Form className="p-3">
            <LoadingOverlay active={loading} />
            <Row>
              <Col lg="3">
                <Select
                  id="kataegisSamples"
                  label="Sample Name"
                  value={sample}
                  options={sampleOptions}
                  onChange={(sample) =>
                    dispatchKataegis({
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
                      dispatchKataegis({ highlight: !highlight });
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
                      dispatchKataegis({
                        min: e.target.value,
                      });
                    }}
                  />
                </Group>
              </Col>
              <Col lg="2">
                <Group controlId="kataegisMax">
                  <Label>Maximum Distance</Label>
                  <Control
                    value={max}
                    placeholder=""
                    onChange={(e) => {
                      dispatchKataegis({
                        max: e.target.value,
                      });
                    }}
                  />
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
                    dispatchKataegis({
                      chromosome: chromosome,
                    })
                  }
                />
              </Col>
              <Col lg="2" className="d-flex justify-content-end">
                <Button
                  className="mt-auto mb-3"
                  variant="primary"
                  onClick={calculateR}
                >
                  Calculate
                </Button>
              </Col>
            </Row>
          </Form>

          <div id="kataegisPlot">
            {err && (
              <div>
                <p>Error: {err}</p>
              </div>
            )}
            <div style={{ display: plotURL ? 'block' : 'none' }}>
              <hr />
              <Plot
                className="p-3"
                plotName={plotPath.split('/').slice(-1)[0]}
                plotURL={plotURL}
                txtPath={txtPath ? projectID + txtPath : null}
              />
              <KataegisTable />
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
