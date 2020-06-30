import React from 'react';
import { useSelector } from 'react-redux';
import { Modal, Button } from 'react-bootstrap';
import { store, updateError } from '../../../services/store';

export function ErrorModal(props) {
  const error = useSelector((store) => store.error);
  // const dispatch = useDispatch();
  const closeErrorModal = () => store.dispatch(updateError({param: 'visible', data: false}));

  return (
    <Modal
      data-testid="ErrorModal"
      show={error.visible}
      onHide={closeErrorModal}
    >
      <Modal.Header closeButton>
        <Modal.Title>Error</Modal.Title>
      </Modal.Header>

      <Modal.Body>
        <p
          data-testid="ErrorModalMessage"
          dangerouslySetInnerHTML={{ __html: error.message }}
        />
      </Modal.Body>

      <Modal.Footer>
        <Button variant="secondary" onClick={closeErrorModal}>
          Close
        </Button>
      </Modal.Footer>
    </Modal>
  );
}
