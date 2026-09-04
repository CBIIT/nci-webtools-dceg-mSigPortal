import { Suspense, useCallback } from 'react';
import { Alert, Button, Container } from 'react-bootstrap';
import { useDispatch, useSelector } from 'react-redux';
import { useRecoilValue } from 'recoil';
import { actions } from '../../../../services/store/visualization';
import Loader from '../../../controls/loader/loader';
import ErrorBoundary from '../../../controls/errorBoundary/error-boundary';
import D3TreeLeaf from './treeLeafPlot';
import TreeLeafForm from './treeLeafForm';
import { exportSvg } from './treeLeaf.utils';
import { formState } from './treeLeaf.state';

export default function TreeAndLeaf({ state = {}, ...props }) {
  const dispatch = useDispatch();
  const publicForm = useSelector((store) => store.visualization.publicForm);
  const { source } = state;
  const form = useRecoilValue(formState);
  const plotId = 'treeLeafPlot';
  const isUser = source === 'user';

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
            sample: event.Sample,
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
          <div className="d-flex justify-content-between align-items-end">
            <TreeLeafForm state={state} />
            <Button variant="link" onClick={handleExport}>
              Export Plot
            </Button>
          </div>
          <D3TreeLeaf
            id={plotId}
            width={2000}
            height={2000}
            onSelect={handleSelect}
            state={state}
          />
        </Suspense>
      </ErrorBoundary>
    </Container>
  );
}
