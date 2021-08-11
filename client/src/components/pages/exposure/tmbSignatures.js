import React from 'react';
import { useSelector } from 'react-redux';
import Description from '../../controls/description/description';
import Plot from '../../controls/plot/plot';
import Debug from '../../controls/debug/debug';

export default function TmbSignatures() {
  const exposure = useSelector((state) => state.exposure);
  const { plotPath, debugR, err } = exposure.tmbSignatures;
  const { projectID } = exposure.exposureState;

  return (
    <div>
      <div className="p-3">
        <b>Tumor Mutational Burden Separated by Signatures</b>
        <Description
          className="m-0"
          short={
            'The bar plot below illustrates the level of tumor mutational burden across different signatures from the selected signature set in the selected cancer type.'
          }
          description={
            'The bar plot below illustrates the level of tumor mutational burden across different signatures from the selected signature set in the selected cancer type. Across the top of the plot are the signatures (in the selected cancer type) that were exhibited by at least one sample. On the y-axis is the number of mutations per Megabase (log10), and the x-axis denotes sample numbers. The green number is the total number of samples with the selected cancer type (and therefore evaluated for each signature). The blue number is the number of samples with the selected cancer type that detected the mutational signature from the reference signature set.'
          }
        />{' '}
      </div>
      <hr />
      <div id="tmbSigPlot">
        {err && (
          <div>
            <hr />
            <p className="p-3 text-danger">{err}</p>
          </div>
        )}
        {plotPath && (
          <>
            <hr />

            <Plot
              className="p-3"
              title="Tumor Mutational Burden Separated by Signatures"
              downloadName={plotPath.split('/').slice(-1)[0]}
              plotPath={`api/results/${projectID}${plotPath}`}
            />
          </>
        )}
      </div>
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}