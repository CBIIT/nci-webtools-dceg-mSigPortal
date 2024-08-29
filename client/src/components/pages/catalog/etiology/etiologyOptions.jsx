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

  // custom sort rules for etiology options
  // order "CN:" prefixed etiologies to the end
  // custom order for "Unknown" etiology
  // otherwise use natural sort
  function etiologySort(a, b) {
    const cnRegex = /^CN:/;
    const unkRegex = /^Unknown/;
    const unkOrder = [
      'Unknown chemotherapy treatment',
      'Unknown (clock-like signature)',
      'Unknown',
    ];

    if (a.match(cnRegex)) {
      return 1;
    } else if (b.match(cnRegex)) {
      return -1;
    } else if (a.match(unkRegex) || b.match(unkRegex)) {
      return unkOrder.indexOf(a) - unkOrder.indexOf(b);
    } else if (a)
      return a.localeCompare(b, undefined, {
        numeric: true,
        sensitivity: 'base',
      });
    else {
      return 0;
    }
  }

  return (
    <div className="mb-3">
      <h5 className="separator">{data[0].json.etiologyDisplay}</h5>
      {data ? (
        <Row className="justify-content-center">
          {[...new Set(data.map((e) => e.etiology))]
            .sort(etiologySort)
            .map((e) => (
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
