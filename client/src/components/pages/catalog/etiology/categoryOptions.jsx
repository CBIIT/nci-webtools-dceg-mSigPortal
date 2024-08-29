import { Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../../services/store/catalog';

import './etiology.scss';

const actions = { ...catalogActions };

export default function CategoryOptions({ categories }) {
  const dispatch = useDispatch();
  const mergeEtiology = (state) =>
    dispatch(actions.mergeCatalog({ etiology: state }));

  const { category } = useSelector((state) => state.catalog.etiology);

  function handleCatalogClick(e) {
    if (e.category != category) {
      mergeEtiology({
        category: e.category,
        etiology: '',
        signature: '',
        study: '',
      });
    }
  }

  return (
    <Row className="justify-content-center mb-3">
      {categories.map((e) => (
        <Col key={e.name} lg="3" md="3" sm="4" className="mb-3 d-flex">
          <Button
            size="sm"
            variant="dark"
            onClick={() => handleCatalogClick(e)}
            className={category != e.category ? 'disabled' : ''}
            block
          >
            <div>
              <div>{e.name}</div>
              <div>{e.author}</div>
            </div>
          </Button>
        </Col>
      ))}
    </Row>
  );
}
