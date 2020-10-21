import React from 'react';
import { useSelector } from 'react-redux';
import { Modal, Button } from 'react-bootstrap';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faTimesCircle } from '@fortawesome/free-solid-svg-icons';
import { store, updateError } from '../../../services/store';

export function ErrorModal(props) {
  const error = useSelector((store) => store.error);
  const closeErrorModal = () => store.dispatch(updateError({ visible: false }));

  return (
    <Modal
      data-testid="ErrorModal"
      show={error.visible}
      onHide={closeErrorModal}
    >
      <Modal.Header style={{ backgroundColor: '#d32f2f' }}>
        <Modal.Title className="d-flex justify-content-center w-100">
          <FontAwesomeIcon
            icon={faTimesCircle}
            size="3x"
            style={{ color: '#fafafa' }}
          />
        </Modal.Title>
      </Modal.Header>

      <Modal.Body
        className="d-flex justify-content-center"
        style={{ backgroundColor: '#fafafa' }}
      >
        <p
          className="m-0 w-100"
          data-testid="ErrorModalMessage"
          dangerouslySetInnerHTML={{ __html: error.message }}
        />
      </Modal.Body>

      <Modal.Footer
        className="d-flex justify-content-center border-0"
        style={{ backgroundColor: '#fafafa' }}
      >
        <Button variant="outline-secondary" onClick={closeErrorModal}>
          Close
        </Button>
      </Modal.Footer>
    </Modal>
  );
}
