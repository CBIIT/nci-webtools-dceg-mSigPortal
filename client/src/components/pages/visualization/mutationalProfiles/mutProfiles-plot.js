import { useEffect, useState } from 'react';
import cloneDeep from 'lodash/cloneDeep';
import { Button, Row, Col } from 'react-bootstrap';
import Plot from 'react-plotly.js';
import { downloadImage } from 'plotly.js';
import { saveAs } from 'file-saver';
import { useSelector } from 'react-redux';
import { useSignatureQuery } from '../../../../services/apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import SBS96 from '../../../controls/plotly/mutationalSignature/sbs96';
import SBS192 from '../../../controls/plotly/mutationalSignature/sbs192';
import SBS1536 from '../../../controls/plotly/mutationalSignature/sbs1536';
import DBS78 from '../../../controls/plotly/mutationalSignature/dbs78';
import ID83 from '../../../controls/plotly/mutationalSignature/id83';

export default function MutProfilePlot() {
  const store = useSelector((state) => state.visualization);

  const { sample, profile, matrix } = store.mutationalProfiles;
  const { study, cancer, strategy } = store.publicForm;

  const [params, setParams] = useState(null);
  const [plot, setPlot] = useState(null);

  const {
    data = [],
    error,
    isFetching,
  } = useSignatureQuery(params, {
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
      };
      setParams(params);
    }
  }, [sample, profile, matrix]);

  // generate plot
  useEffect(() => {
    if (data.length) generatePlot(data);
  }, [data]);

  function generatePlot(data) {
    const profileMatrix = profile.value + matrix.value;

    const { traces, layout } =
      profileMatrix == "SBS96"
        ? SBS96(data)
        : profileMatrix == "SBS192"
        ? SBS192(data)
        : profileMatrix == "SBS1536"
        ? SBS1536(data)
        : profileMatrix == "DBS78"
        ? DBS78(data)
        : profileMatrix == "ID83"
        ? ID83(data)
        : { traces: [], layout: {} };

    setPlot({ data: [...traces], layout });
  }

  const divId = "mutationalProfilePlot";
  const config = {
    displayModeBar: true,
    responsive: true,
    displaylogo: false,
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
            <Plot
              className="w-100"
              divId={divId}
              style={{ height: "500px" }}
              data={cloneDeep(plot.data)}
              layout={cloneDeep(plot.layout)}
              config={cloneDeep(config)}
              useResizeHandler
            />
            <Row className="justify-content-center">
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
            </Row>
          </div>
        )
      )}
    </>
  );
}
