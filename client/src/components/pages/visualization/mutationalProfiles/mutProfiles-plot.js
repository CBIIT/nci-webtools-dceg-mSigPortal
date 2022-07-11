import { useEffect, useState } from 'react';
import cloneDeep from 'lodash/cloneDeep';
import { Button, Row, Col } from 'react-bootstrap';
import Plot from 'react-plotly.js';
import { downloadImage } from 'plotly.js';
import { saveAs } from 'file-saver';
import { useSelector } from 'react-redux';
import { useSeqmatrixQuery } from '../../../../services/apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import SBS6 from '../../../controls/plotly/mutationalSignature/sbs6';
import SBS24 from '../../../controls/plotly/mutationalSignature/sbs24';
import SBS96 from '../../../controls/plotly/mutationalSignature/sbs96';
import SBS192 from '../../../controls/plotly/mutationalSignature/sbs192';
import SBS288 from '../../../controls/plotly/mutationalSignature/sbs288';
import SBS384 from '../../../controls/plotly/mutationalSignature/sbs384';
import SBS1536 from '../../../controls/plotly/mutationalSignature/sbs1536';
import DBS78 from '../../../controls/plotly/mutationalSignature/dbs78';
import DBS186 from '../../../controls/plotly/mutationalSignature/dbs186';
import ID83 from '../../../controls/plotly/mutationalSignature/id83';
import ID28 from '../../../controls/plotly/mutationalSignature/id28';
import ID415 from '../../../controls/plotly/mutationalSignature/id415';

export default function MutProfilePlot() {
  const store = useSelector((state) => state.visualization);

  const { sample, profile, matrix } = store.mutationalProfiles;
  const { study, cancer, strategy } = store.publicForm;
  const { source } = store.main;

  const [params, setParams] = useState(null);
  const [plot, setPlot] = useState(null);

  const {
    data = [],
    error,
    isFetching,
  } = useSeqmatrixQuery(params, {
    skip: !params,
  });

  // get data on form change
  useEffect(() => {
    if (sample && sample.value && profile.value && matrix.value) {
      const params = {
        study: study.value,
        cancer: cancer.value,
        strategy: strategy.value,
        sample: sample.value,
        profile: profile.value,
        matrix: matrix.value,
        type: 'mutationalProfiles',
      };
      setParams(params);
    }
  }, [sample, profile, matrix]);

  // generate plot
  useEffect(() => {
    if (data.length) generatePlot(data, sample.value);
  }, [data]);

  function generatePlot(data, sample) {
    const profileMatrix = profile.value + matrix.value;
    const { traces, layout } =
      profileMatrix == 'SBS6'
        ? SBS6(data, sample)
        : profileMatrix == 'SBS24'
        ? SBS24(data, sample)
        : profileMatrix == 'SBS96'
        ? SBS96(data, sample)
        : profileMatrix == 'SBS192'
        ? SBS192(data, sample)
        : profileMatrix == 'SBS288'
        ? SBS288(data, sample)
        : profileMatrix == 'SBS384'
        ? SBS384(data, sample)
        : profileMatrix == 'SBS1536'
        ? SBS1536(data, sample)
        : profileMatrix == 'DBS78'
        ? DBS78(data, sample)
        : profileMatrix == 'DBS186'
        ? DBS186(data, sample)
        : profileMatrix == 'ID28'
        ? ID28(data, sample)
        : profileMatrix == 'ID83'
        ? ID83(data, sample)
        : profileMatrix == 'ID415'
        ? ID415(data, sample)
        : { traces: [], layout: {} };

    setPlot({ data: [...traces], layout });
  }

  const divId = 'mutationalProfilePlot';
  const config = {
    displayModeBar: true,
    responsive: true,
    displaylogo: false,
    toImageButtonOptions: {
      format: 'svg',
      filename: 'plot_export',
      height: 1000,
      width: 1000,
      scale: 1,
    },
  };

  return (
    <>
      <LoadingOverlay active={isFetching} />
      {error ? (
        <div className="text-center">
          <div>An error as occured</div>
          <div>{error.data}</div>
        </div>
      ) : (
        plot && (
          <div className="mb-3">
            <Row className="justify-content-center">
              <Plot
                {...(profile.value + matrix.value === 'SBS1536' ||
                profile.value + matrix.value === 'SBS6' ||
                profile.value + matrix.value === 'SBS24' ||
                profile.value + matrix.value === 'ID28'
                  ? { className: 'w-70' }
                  : { className: 'w-95' })}
                divId={divId}
                style={{
                  height: '600px',
                }}
                data={cloneDeep(plot.data)}
                layout={cloneDeep(plot.layout)}
                config={cloneDeep(config)}
                useResizeHandler
              />
            </Row>

            {/* <Row className="justify-content-center">
              <Col sm="auto">
                <Button
                  onClick={() =>
                    downloadImage(divId, {
                      format: "png",
                      filename: sample.value,
                    })
                  }
                >
                  Download PNG
                </Button>
              </Col>
              <Col sm="auto">
                <Button
                  onClick={() =>
                    downloadImage(divId, {
                      format: "svg",
                      filename: sample.value,
                    })
                  }
                >
                  Download SVG
                </Button>
              </Col>
              <Col sm="auto">
                <Button
                  onClick={() =>
                    saveAs(
                      new Blob([JSON.stringify(plot)], {
                        type: "application/json",
                      }),
                      `${sample.value}.json`
                    )
                  }
                >
                  Download JSON
                </Button>
              </Col>
            </Row> */}
          </div>
        )
      )}
    </>
  );
}
