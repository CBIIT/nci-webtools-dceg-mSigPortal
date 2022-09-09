import { Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../../services/store/catalog';
import { actions as modalActions } from '../../../../services/store/modal';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

const actions = { ...catalogActions, ...modalActions };

export default function EtiologyOptions({ data }) {
  const dispatch = useDispatch();
  const mergeEtiology = (state) =>
    dispatch(actions.mergeCatalog({ etiology: state }));

  const { category, etiology } = useSelector((state) => state.catalog.etiology);

  return (
    <div className="mb-3">
      <h5 className="separator">{data[0].etiologyDisplay}</h5>
      {data ? (
        <Row className="justify-content-center">
          {[...new Set(data.map((e) => e.etiology))].sort().map((e) => (
            <Col
              key={e}
              lg={category == 'GeneEdits' ? '1' : '2'}
              md="3"
              sm="4"
              className="mb-3 d-flex"
            >
              <Button
                size="sm"
                variant="dark"
                onClick={() =>
                  mergeEtiology({
                    etiology: e,
                    signature: '',
                  })
                }
                className={etiology != e ? 'disabled' : ''}
                block
              >
                {e}
              </Button>
            </Col>
          ))}
        </Row>
      ) : (
        <div style={{ minHeight: '100px' }}>
          <LoadingOverlay active={true} />
        </div>
      )}
    </div>
  );
}
