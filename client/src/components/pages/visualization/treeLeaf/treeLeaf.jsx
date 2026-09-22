import { Suspense, useCallback, useState } from 'react';
import { Alert, Button, Container } from 'react-bootstrap';
import { useDispatch, useSelector } from 'react-redux';
import { useRecoilValue } from 'recoil';
import { actions } from '@/services/store/visualization';
import Loader from '@/components/controls/loader/loader';
import ErrorBoundary from '@/components/controls/errorBoundary/error-boundary';
import D3TreeLeaf from './treeLeafPlot';
import TreeLeafForm from './treeLeafForm';
import { exportSvg } from './treeLeaf.utils';
import {
  graphDataSelector,
  getInitialFormState,
  SIGNATURE_SET_NAME,
  PROFILE,
  MATRIX,
} from './treeLeaf.state';

const plotId = 'treeLeafPlot';

export default function TreeAndLeaf({ state = {}, ...props }) {
  const publicForm = useSelector((store) => store.visualization.publicForm);
  const isUser = state.source === 'user';

  const resetKey = isUser
    ? `user:${state.id ?? ''}`
    : `public:${publicForm?.study?.value ?? ''}:${
        publicForm?.cancer?.value ?? ''
      }`;

  return (
    <TreeLeafView key={resetKey} state={state} isUser={isUser} {...props} />
  );
}

function TreeLeafView({ state, isUser, ...props }) {
  const dispatch = useDispatch();
  const publicForm = useSelector((store) => store.visualization.publicForm);
  const [form, setForm] = useState(() =>
    getInitialFormState({ isUser, publicForm })
  );

  const mergeForm = useCallback(
    (next) => setForm((prev) => ({ ...prev, ...next })),
    []
  );

  function handleExport() {
    const plotSelector = `#${plotId}`;
    const studyLabel = isUser ? 'User Data' : publicForm?.study?.label;
    const fileName = `treeLeafPlot ${studyLabel} ${form.color.label}.svg`;
    exportSvg(plotSelector, fileName);
  }

  const handleSelect = useCallback(
    (event) => {
      dispatch(
        actions.mergeVisualization({
          main: {
            displayTab: 'mutationalProfiles',
            openSidebar: false,
          },
          mutationalProfiles: {
            sample: event.SampleName ?? event.Sample,
            filter: event.Filter ?? '',
          },
        })
      );
    },
    [dispatch]
  );

  const defaultFallback = isUser
    ? 'The selected user session does not provide SBS96 mutation seqmatrix data.'
    : 'The selected study does not provide exposure and mutation seqmatrix data.';

  return (
    <Container
      fluid
      className="bg-white border rounded p-3 text-center"
      style={{ minHeight: 500 }}
      {...props}
    >
      <ErrorBoundary
        fallback={(error) => (
          <Alert variant="warning">{error?.message || defaultFallback}</Alert>
        )}
      >
        <Suspense fallback={<Loader message="Loading Plot Data" />}>
          <TreeLeafContent
            state={state}
            isUser={isUser}
            form={form}
            onChange={mergeForm}
            onExport={handleExport}
            onSelect={handleSelect}
          />
        </Suspense>
      </ErrorBoundary>
    </Container>
  );
}

function TreeLeafContent({
  state,
  isUser,
  form,
  onChange,
  onExport,
  onSelect,
}) {
  const publicForm = useSelector((store) => store.visualization.publicForm);
  const { id: sessionId } = state;
  const study = publicForm?.study?.value ?? 'PCAWG';
  const strategy = publicForm?.strategy?.value ?? 'WGS';

  const params = isUser
    ? { userId: sessionId, source: 'user', profile: PROFILE, matrix: MATRIX }
    : {
        study,
        strategy,
        signatureSetName: SIGNATURE_SET_NAME,
        profile: PROFILE,
        matrix: MATRIX,
        cancer: form.cancer,
      };

  const graphData = useRecoilValue(graphDataSelector(params));
  if (graphData?.error) {
    throw new Error(graphData.error);
  }
  // Trigger ErrorBoundary when the request fails or returns no data
  if (graphData === null && !(isUser && !sessionId)) {
    throw new Error(
      isUser
        ? 'Failed to load Tree and Leaf data for this user session.'
        : 'Failed to load Tree and Leaf data for the selected study.'
    );
  }

  const { attributes } = graphData || {};

  return (
    <>
      <div className="d-flex justify-content-between align-items-end">
        <TreeLeafForm
          isUser={isUser}
          form={form}
          onChange={onChange}
          attributes={attributes}
        />
        <Button variant="link" onClick={onExport}>
          Export Plot
        </Button>
      </div>
      <D3TreeLeaf
        id={plotId}
        width={2000}
        height={2000}
        onSelect={onSelect}
        isUser={isUser}
        graphData={graphData}
        form={form}
      />
    </>
  );
}
