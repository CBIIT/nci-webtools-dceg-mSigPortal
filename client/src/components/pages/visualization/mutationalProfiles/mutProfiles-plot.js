import cloneDeep from "lodash/cloneDeep";
import { Button, Row, Col } from "react-bootstrap";
import Plot from "react-plotly.js";
import { downloadImage } from "plotly.js";
import { saveAs } from "file-saver";
import { useSelector, useDispatch } from "react-redux";
import { actions as visualizationActions } from "../../../../services/store/visualization";
import SBS96 from "../../../controls/plotly/mutationalSignature/sbs96";
import SBS192 from "../../../controls/plotly/mutationalSignature/sbs192";
import { useEffect } from "react";
import axios from "axios";
import DBS78 from "../../../controls/plotly/mutationalSignature/dbs78";
import ID83 from "../../../controls/plotly/mutationalSignature/id83";
import { LoadingOverlay } from "../../../controls/loading-overlay/loading-overlay";

export default function MutProfilePlot() {
  const dispatch = useDispatch();
  const store = useSelector((state) => state.visualization);

  const mergeState = (state) =>
    dispatch(
      visualizationActions.mergeVisualization({ mutationalProfiles: state })
    );

  const { sample, profile, matrix, data, plot, loading } =
    store.mutationalProfiles;
  const { study, cancer, strategy } = store.publicForm;

  // get data
  useEffect(() => {
    if (sample && sample.value && profile.value && matrix.value)
      getSignatureData();
  }, [sample, profile, matrix]);

  // generate plot
  useEffect(() => {
    if (data.length) generatePlot(data);
  }, [data]);

  async function getSignatureData() {
    mergeState({ loading: true });
    try {
      const { data } = await axios.get("web/signature", {
        params: {
          study: study.value,
          cancer: cancer.value,
          strategy: strategy.value,
          sample: sample.value,
          profile: profile.value,
          matrix: matrix.value,
        },
      });
      mergeState({ data: data });
    } catch (error) {}
    mergeState({ loading: false });
  }

  function generatePlot(data) {
    mergeState({ loading: true });
    const profileMatrix = profile.value + matrix.value;
    console.log(profileMatrix);
    const { traces, layout } =
      profileMatrix == "SBS96"
        ? SBS96(data)
        : profileMatrix == "SBS192"
        ? SBS192(data)
        : profileMatrix == "DBS78"
        ? DBS78(data)
        : profileMatrix == "ID83"
        ? ID83(data)
        : { traces: [], layout: {} };

    mergeState({
      plot: { data: [...traces], layout },
      loading: false,
    });
  }

  const divId = "mutationalProfilePlot";
  const config = {
    displayModeBar: true,
    response: true,
    displaylogo: false,
  };

  return (
    <div>
      <LoadingOverlay active={loading} />
      {plot && plot.data.length ? (
        <div>
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
      ) : (
        <div>Select a supported profile</div>
      )}
    </div>
  );
}
