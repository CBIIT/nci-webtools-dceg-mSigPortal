import { Suspense, useCallback, useMemo, useState } from 'react';
import { Alert, Button, Container } from 'react-bootstrap';
import { useDispatch, useSelector } from 'react-redux';
import { useRecoilValue } from 'recoil';
import { actions } from '@/services/store/visualization';
import Loader from '@/components/controls/loader/loader';
import ErrorBoundary from '@/components/controls/errorBoundary/error-boundary';
import D3TreeLeaf from './treeLeafPlot';
import TreeLeafForm from './treeLeafForm';
import { exportSvg, groupBy, getPlotTitle } from './treeLeaf.utils';
import {
  graphDataSelector,
  getInitialFormState,
  SIGNATURE_SET_NAME,
  PROFILE,
  MATRIX,
} from './treeLeaf.state';

const plotId = 'treeLeafPlot';
// the server lays out the tree with radius plotSize / 2 (TREE_LEAF_RADIUS)
const plotSize = 2000;

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
  const { nodes, links, attributes, params: parameters } = graphData || {};
  const attributesBySample = useMemo(
    () => (attributes ? groupBy(attributes, 'Sample') : null),
    [attributes]
  );
  // nodes/links arrive already laid out by the server; just pass them through
  const layout = useMemo(
    () => (nodes && links ? { nodes, links } : null),
    [nodes, links]
  );

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

  const plotTitle = getPlotTitle({
    isUser,
    form,
    studyLabel: publicForm?.study?.label,
    signatureSetName: parameters?.signatureSetName,
  });

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
      <div className="border rounded p-3 position-relative">
        <D3TreeLeaf
          id={plotId}
          width={plotSize}
          height={plotSize}
          onSelect={onSelect}
          layout={parameters ? layout : null}
          attributes={attributesBySample}
          form={form}
          plotTitle={plotTitle}
        />
      </div>
    </>
  );
}
