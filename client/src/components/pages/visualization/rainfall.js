import React, { useEffect } from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { dispatchError, dispatchRainfall } from '../../../services/store';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import Plot from '../../controls/plot/plot';
import Debug from '../../controls/debug/debug';
import Select from '../../controls/select/select';
import Accordions from '../../controls/accordions/accordions';

const { Group, Check, Label, Control } = Form;

export default function PCA({ submitR }) {
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
  } = useSelector((state) => state.rainfall);

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
          dispatchRainfall({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchRainfall({ err: true, plotURL: '' });
    }
  }

  function clearPlot() {
    if (plotURL) {
      URL.revokeObjectURL(plotURL);
      dispatchRainfall({ plotPath: '', plotURL: '' });
    }
  }

  async function calculateR() {
    dispatchRainfall({
      loading: true,
      err: false,
      debugR: '',
    });

    try {
      const response = await submitR('rainfall', {
        sample: sample,
        highlight: highlight,
        min: parseInt(min),
        max: parseInt(max),
        chromosome: chromosome,
      });
      if (!response.ok) {
        const err = await response.json();
        dispatchRainfall({ debugR: err });

        dispatchRainfall({ loading: false });
      } else {
        const { debugR, output, errors } = await response.json();
        if (Object.keys(output).length) {
          dispatchRainfall({
            debugR: debugR,
            plotPath: output.plotPath,
            // txtPath: output.txtPath,
            loading: false,
          });
        } else {
          dispatchRainfall({
            debugR: debugR,
            plotPath: '',
            // txtPath: '',
            err: errors,
            loading: false,
          });
        }
      }
    } catch (err) {
      dispatchError(err);
      dispatchRainfall({ loading: false });
    }
  }

  const plots = (
    <div id="rainfallPlot">
      {err && (
        <div>
          <p>Error: {err}</p>
        </div>
      )}
      <div style={{ display: plotURL ? 'block' : 'none' }}>
        <Plot
          plotName={plotPath.split('/').slice(-1)[0]}
          plotURL={plotURL}
          txtPath={projectID + txtPath}
        />
      </div>
    </div>
  );

  const component =
    source == 'user' && inputFormat == 'vcf' ? (
      <Form>
        <LoadingOverlay active={loading} />
        <div>
          <Row className="">
            <Col lg="3">
              <Select
                id="rainfallSamples"
                label="Sample Name"
                value={sample}
                options={sampleOptions}
                onChange={(sample) =>
                  dispatchRainfall({
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
                    dispatchRainfall({ highlight: !highlight });
                  }}
                />
              </Group>
            </Col>
            <Col lg="2">
              <Group controlId="rainfallMin">
                <Label>Minimum Number of Mutations</Label>
                <Control
                  value={min}
                  placeholder=""
                  onChange={(e) => {
                    dispatchRainfall({
                      min: e.target.value,
                    });
                  }}
                />
                {/* <Text className="text-muted">(Ex. NCG>NTG)</Text> */}
              </Group>{' '}
            </Col>
            <Col lg="2">
              <Group controlId="rainfallMax">
                <Label>Maximum Distance</Label>
                <Control
                  value={max}
                  placeholder=""
                  onChange={(e) => {
                    dispatchRainfall({
                      max: e.target.value,
                    });
                  }}
                />
                {/* <Text className="text-muted">(Ex. NCG>NTG)</Text> */}
              </Group>
            </Col>
            <Col lg="2">
              <Select
                id="rainfallChromosome"
                label="Chromosome"
                value={chromosome}
                options={['None', ...[...Array(23).keys()].slice(1), 'X', 'Y']}
                onChange={(chromosome) =>
                  dispatchRainfall({
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
          {plots}
        </div>
      </Form>
    ) : (
      <p>
        Kataegis Identification is only available for Data Source: User and
        requires VCF input files
      </p>
    );

  const accordions = [
    { title: 'Kataegis Identification', component: component },
  ];

  return (
    <div>
      <Accordions components={accordions} />
      <Debug msg={debugR} />
    </div>
  );
}
