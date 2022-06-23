import React from 'react';
import { useSelector } from 'react-redux';
import Description from '../../controls/description/description';
import SvgContainer from '../../controls/svgContainer/svgContainer';
import Debug from '../../controls/debug/debug';

export default function TMB() {
  const exposure = useSelector((state) => state.exposure);
  const { plotPath, debugR, err } = exposure.tmb;
  const { projectID } = exposure.exposureState;

  return (
    <div>
      <div className="p-3">
        <b>Tumor Mutational Burden</b>
        <Description
          className="m-0"
          less="The bar plot below illustrates tumor mutational burden (number of mutations per Mb) across different cancer types for the selected study."
          more="The y-axis is the number of mutations per Mb (log10), and the x-axis indicates cancer types. The green number is the number of samples for a given cancer type, and the blue number is the number of samples that have mutation data available for that cancer type."
        />
      </div>
      <div id="tmbPlot">
        {err && (
          <div>
            <hr />
            <p className="p-3 text-danger">{err}</p>
          </div>
        )}
        {plotPath && (
          <>
            <hr />
            <SvgContainer
              className="p-3"
              title="Tumor Mutational Burden"
              downloadName={plotPath.split('/').slice(-1)[0]}
              plotPath={`web/results/${plotPath}`}
            />
          </>
        )}
      </div>

      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
