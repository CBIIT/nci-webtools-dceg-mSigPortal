import React from 'react';
import { useSelector } from 'react-redux';
import Description from '../../controls/description/description';
import Plot from '../../controls/plot/plot';

export default function MsDecomposition() {
  const exposure = useSelector((state) => state.exposure);
  const { plotPath, txtPath, debugR, err } = exposure.msDecomposition;

  return (
    <div>
      <div className="p-3">
        <b>Evaluating the Performance of Mutational Signature Decomposition</b>
        <Description
          less="The distribution plot below illustrates mutational signature decomposition distribution in a selected cancer type (by selecting [Cancer Type Only] on the left panel) or across different cancer types."
          more="Five different methods are used to measure the similarities of original mutational profile matrix and reconstructed mutational profile matrix across all the samples: cosine similarity, 100-L1_Norm_%, 100-L2_Norm_%, KL_Divergence, and Pearson Correlation. Click here for details of these evaluation methods."
        />
      </div>
      <hr />

      <div id="decompositionPlot">
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
              title="Evaluating the Performance of Mutational Signature Decomposition"
              downloadName={plotPath.split('/').slice(-1)[0]}
              plotPath={`api/results/${plotPath}`}
              txtPath={`api/results/${txtPath}`}
            />
          </>
        )}
      </div>
      {/* <Debug msg={debugR} /> */}
    </div>
  );
}
