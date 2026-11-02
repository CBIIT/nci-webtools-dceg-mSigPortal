import { useState, useEffect } from 'react';
import Plotly from '@/components/controls/plotly/plot/plot';
import { useMsBurdenQuery } from './apiSlice';
import { LoadingOverlay } from '@/components/controls/loading-overlay/loading-overlay';
import { getErrorMessage, getDatasetLabel } from '@/services/utils';

export default function MsBurdenPlot({
  state,
  form,
  availableOptions,
  divId = 'msBurden',
}) {
  const [params, setParams] = useState('');
  const { currentData, error, isFetching } = useMsBurdenQuery(params, {
    skip: !params,
  });

  const { signatureName } = form;
  const { study, strategy, signatureSetName, cancer, useAllCancer, id } = state;

  const signatureMissing =
    signatureName &&
    availableOptions &&
    !availableOptions.some((e) => e.value === signatureName.value);

  // query after signature name is changed
  useEffect(() => {
    if (signatureMissing) return;
    if (signatureName && id) {
      setParams({
        signatureName: signatureName.value,
        userId: id,
      });
    } else if (signatureName && study) {
      setParams({
        signatureName: signatureName.value,
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        ...(!useAllCancer && { cancer: cancer.value }),
      });
    }
  }, [signatureName, id, signatureMissing]);

  const label = getDatasetLabel(state);

  return (
    <div>
      <LoadingOverlay active={isFetching} />
      <h5 className="d-flex justify-content-center">
        Mutational Signature Burden Across Cancer Types
      </h5>
      {label && <p className="text-center text-muted mb-0">{label}</p>}
      {signatureMissing ? (
        <p className="p-3 text-danger text-center">
          Signature {signatureName.value} is missing in this dataset
        </p>
      ) : (
        <>
          {error && <p className="p-3 text-danger">{getErrorMessage(error)}</p>}
          {currentData && (
            <Plotly
              className="w-100"
              divId={divId}
              filename={divId}
              data={currentData.traces}
              layout={currentData.layout}
              config={currentData.config}
            />
          )}
        </>
      )}
    </div>
  );
}
