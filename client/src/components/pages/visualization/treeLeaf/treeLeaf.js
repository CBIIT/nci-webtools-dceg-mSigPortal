import { Suspense } from 'react';
import { Alert, Button, Container } from 'react-bootstrap';
import Loader from '../../../controls/loader/loader';
import ErrorBoundary from '../../../controls/errorBoundary/error-boundary';
import D3TreeLeaf from './treeLeafPlot';
import TreeLeafForm from './treeLeafForm';
import { exportSvg } from './treeLeaf.utils';

export default function TreeAndLeaf(props) {
  const plotId = 'treeLeafPlot';
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
            <Button variant="link" onClick={() => exportSvg(`#${plotId}`, 'treeLeaf.svg')}>Export Plot</Button>
          </div>
          <D3TreeLeaf id={plotId} width={2000} height={2000} onSelect={props.onSelect} />
        </Suspense>
      </ErrorBoundary>
    </Container>
  );
}
