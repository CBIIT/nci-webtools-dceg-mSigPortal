
import { Suspense } from 'react';
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


export default function TreeAndLeaf(props) {
  const dispatch = useDispatch();
  const mergeVisualizationState = (state) => dispatch(actions.mergeVisualization({ main: state }));
  const publicForm = useSelector((state) => state.visualization.publicForm);
  const form = useRecoilValue(formState);
  const plotId = 'treeLeafPlot';

  function handleExport() {
    const plotSelector = `#${plotId}`;
    const fileName = `treeLeafPlot ${publicForm?.study?.label} ${form.color.label}.svg`;
    exportSvg(plotSelector, fileName);
  }

  function handleSelect(event) {
    dispatch(actions.mergeVisualization({
      main: {
        displayTab: 'mutationalProfiles', 
        openSidebar: false
      },
      mutationalProfiles: {
        sample: event.Sample,
      }
    }))
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
          <Alert variant="warning">
            The selected study does not provide exposure and mutation seqmatrix data.
          </Alert>
        }
      >
        <Suspense fallback={<Loader message="Loading Plot" />}>
          <div className="d-flex justify-content-between align-items-end">
            <TreeLeafForm />
            <Button variant="link" onClick={handleExport}>Export Plot</Button>
          </div>
          <D3TreeLeaf id={plotId} width={2000} height={2000} onSelect={handleSelect} />
        </Suspense>
      </ErrorBoundary>
    </Container>
  );
}
