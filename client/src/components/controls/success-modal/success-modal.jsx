import React from 'react';
import { useSelector, useDispatch } from 'react-redux';
import { Modal, Button } from 'react-bootstrap';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faCheckCircle } from '@fortawesome/free-solid-svg-icons';
import { actions } from '../../../services/store/modal';

export function SuccessModal(props) {
  const dispatch = useDispatch();
  const mergeSuccess = (state) =>
    dispatch(actions.mergeModal({ success: state }));

  const success = useSelector((state) => state.modal['success']);

  const closeModal = () => mergeSuccess({ visible: false });

  return (
    <Modal
      data-testid="SuccessModal"
      show={success.visible}
      onHide={closeModal}
    >
      <Modal.Header style={{ backgroundColor: '#04c585 ' }}>
        <Modal.Title className="d-flex justify-content-center w-100">
          <FontAwesomeIcon
            icon={faCheckCircle}
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
          className="m-0"
          data-testid="ModalMessage"
          dangerouslySetInnerHTML={{ __html: success.message }}
        />
      </Modal.Body>

      <Modal.Footer
        className="d-flex justify-content-center border-0"
        style={{ backgroundColor: '#fafafa' }}
      >
        <Button variant="outline-secondary" onClick={closeModal}>
          Close
        </Button>
      </Modal.Footer>
    </Modal>
  );
}
