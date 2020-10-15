import React from 'react';
import { Form, Row, Col, Button } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { dispatchError, dispatchExpTumor } from '../../../../services/store';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Plot from '../../../controls/plot/plot';
import Debug from '../../../controls/debug/debug';
import Select from '../../../controls/select/select';

const { Group, Label, Control, Text } = Form;

export default function Tumor({ submitR }) {
  const rootURL = window.location.pathname;
  const {
    study,
    studyOptions,
    strategy,
    strategyOptions,
    refSignatureSet,
    refSignatureSetOptions,
    genomeSize,
    loading,
  } = useSelector((state) => state.expExposure);
  const { plotPath, plotURL, txtPath, debugR, err, displayDebug } = useSelector(
    (state) => state.expTumor
  );
  const { displayTab, publicDataOptions } = useSelector(
    (state) => state.exploring
  );

  async function setRPlot(plotPath) {
    if (plotPath) {
      try {
        const response = await fetch(`${rootURL}results/${plotPath}`);
        if (!response.ok) {
          // console.log(await response.json());
        } else {
          const pic = await response.blob();
          const objectURL = URL.createObjectURL(pic);

          if (plotURL) URL.revokeObjectURL(plotURL);
          dispatchExpTumor({
            plotURL: objectURL,
          });
        }
      } catch (err) {
        dispatchError(err);
      }
    } else {
      if (plotURL) URL.revokeObjectURL(plotURL);
      dispatchExpTumor({ err: true, plotURL: '' });
    }
  }

  return (
    <div>
      <LoadingOverlay active={loading} />
      {err && (
        <p>An error has occured. Check the debug section for more info.</p>
      )}
      {plotURL && (
        <Plot
          plotName={plotPath.split('/').slice(-1)[0]}
          plotURL={plotURL}
          txtPath={txtPath}
        />
      )}
      <Debug msg={debugR} />
    </div>
  );
}
