import cloneDeep from "lodash/cloneDeep";
import { useRecoilValue } from "recoil";
import { getPlot } from "./tmb.state";
import Plot from "react-plotly.js";

export default function MutProfilePlot() {
  const { data, layout, config } = useRecoilValue(getPlot);

  return (
    <div>
      {data.length ? (
        <Plot
          className="w-100"
          style={{ height: "500px" }}
          data={cloneDeep(data)}
          layout={cloneDeep(layout)}
          config={cloneDeep(config)}
          useResizeHandler
        />
      ) : (
        <div>Select a Cancer Type</div>
      )}
    </div>
  );
}
