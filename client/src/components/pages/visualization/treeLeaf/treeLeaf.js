import { Suspense } from 'react';
import { Alert, Button, Container } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import { useRecoilValue } from 'recoil';
import Loader from '../../../controls/loader/loader';
import ErrorBoundary from '../../../controls/errorBoundary/error-boundary';
import D3TreeLeaf from './treeLeafPlot';
import TreeLeafForm from './treeLeafForm';
import { exportSvg } from './treeLeaf.utils';
import { formState } from './treeLeaf.state';

export default function TreeAndLeaf(props) {
  const plotId = 'treeLeafPlot';
  const publicForm = useSelector((state) => state.visualization.publicForm);
  const form = useRecoilValue(formState);

  function handleExport() {
    const plotSelector = `#${plotId}`;
    const fileName = `treeLeafPlot ${publicForm?.study?.label} ${form.color.label}.svg`;
    exportSvg(plotSelector, fileName);
  }

  return (
    <Container
      fluid
      className="bg-white border rounded p-3 text-center"
      style={{ minHeight: 500 }}
      {...props}
    >
      <ErrorBoundary
        fallback={
          <Alert variant="danger">
            An internal error prevented plots from loading. Please contact the
            website administrator if this problem persists.
          </Alert>
        }
      >
        <Suspense fallback={<Loader message="Loading Plot" />}>
          <div className="d-flex justify-content-between align-items-end">
            <TreeLeafForm />
            <Button variant="link" onClick={handleExport}>Export Plot</Button>
          </div>
          <D3TreeLeaf id={plotId} width={2000} height={2000} onSelect={props.onSelect} />
        </Suspense>
      </ErrorBoundary>
    </Container>
  );
}
