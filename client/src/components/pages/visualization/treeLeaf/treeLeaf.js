import { Suspense } from 'react';
import { Alert, Container } from 'react-bootstrap';
import Loader from '../../../controls/loader/loader';
import ErrorBoundary from '../../../controls/errorBoundary/error-boundary';
import D3TreeLeaf from './treeLeafPlot';
import TreeLeafForm from './treeLeafForm';

export default function TreeAndLeaf(props) {
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
          <TreeLeafForm />
          <D3TreeLeaf width={1000} height={1000} />
        </Suspense>
      </ErrorBoundary>
    </Container>
  );
}
