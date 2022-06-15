import cloneDeep from "lodash/cloneDeep";
import { useRecoilValue } from "recoil";
import { getPlot } from "./tmb.state";
import Plot from "react-plotly.js";

export default function MutProfilePlot() {
  const { data, layout, config } = useRecoilValue(getPlot);
  console.log(data);
  return (
    <div>
      <div>
        {data.length ? (
          <div>
            <h1 className="h5">Tumor Mutational Burden</h1>
            <Plot
              className="w-100"
              style={{ height: "500px" }}
              data={cloneDeep(data)}
              layout={cloneDeep(layout)}
              config={cloneDeep(config)}
              useResizeHandler
            />
          </div>
        ) : (
          <div>Select a Cancer Type</div>
        )}
      </div>
    </div>
  );
}
