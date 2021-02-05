import React from 'react';
import { useSelector, useDispatch } from 'react-redux';
import { Modal, Button } from 'react-bootstrap';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faTimesCircle } from '@fortawesome/free-solid-svg-icons';
import { actions } from '../../../services/store/modal';

export function ErrorModal(props) {
  const dispatch = useDispatch();
  const mergeError = (state) => dispatch(actions.mergeModal({ error: state }));

  const error = useSelector((state) => state.modal['error']);
  const closeErrorModal = () => mergeError({ visible: false });

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
