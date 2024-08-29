import {
  Form,
  Row,
  Col,
  Button,
  OverlayTrigger,
  Popover,
} from 'react-bootstrap';
import { Controller } from 'react-hook-form';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faFolderMinus, faInfoCircle } from '@fortawesome/free-solid-svg-icons';

export default function MsLandscapeForm({
  state,
  variableFile,
  formControl,
  setValue,
  resetField,
}) {
  return (
    <Form className="p-3">
      <Row className="">
        <Col lg="auto">
          <Form.Group controlId="landscape">
            <Form.Label>
              Upload Variable Data{' '}
              <OverlayTrigger
                trigger="hover"
                placement="top"
                overlay={
                  <Popover id="upload-variable-info">
                    <Popover.Content>
                      A text file with a header including two columns of data:
                      Samples and Variable Value
                    </Popover.Content>
                  </Popover>
                }
                rootClose
              >
                <FontAwesomeIcon icon={faInfoCircle} className="btn-link" />
              </OverlayTrigger>
            </Form.Label>
            <div className="d-flex">
              <Controller
                name="variableFile"
                control={formControl}
                render={({ field }) => (
                  <Form.File
                    {...field}
                    className="w-100"
                    value={''}
                    id="variableFile"
                    label={variableFile.name || 'Upload here (optional)'}
                    title={variableFile.name || 'Upload here (optional)'}
                    onChange={(e) => {
                      if (e.target.files.length) {
                        setValue('variableFile', e.target.files[0]);
                      }
                    }}
                    custom
                  />
                )}
              />
              {variableFile.size > 0 && (
                <Button
                  className="ml-1"
                  size="sm"
                  title="Remove"
                  variant="danger"
                  onClick={() => resetField('variableFile')}
                >
                  <FontAwesomeIcon icon={faFolderMinus} size="lg" />
                </Button>
              )}
            </div>
          </Form.Group>
        </Col>
      </Row>
    </Form>
  );
}
