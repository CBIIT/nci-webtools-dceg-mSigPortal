import { useState, useEffect } from 'react';
import {
  Form,
  Row,
  Col,
  Button,
  OverlayTrigger,
  Popover,
} from 'react-bootstrap';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faFolderMinus, faInfoCircle } from '@fortawesome/free-solid-svg-icons';
import Plotly from '../../../controls/plotly/plot/plot';
import { useSelector, useDispatch } from 'react-redux';
import { useMsLandscapePlotQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import { useRef } from 'react';

import MsLandscape from '../../../controls/plotly/msLandscape/msLandscape';
import { readFile, parseMatrix } from '../../../controls/utils/utils';
import { actions } from '../../../../services/store/exposure';

const { Label, Group } = Form;
export default function MsLandscapePlot({ state }) {
  // const [params, setParams] = useState('');
  // const { data, error, isFetching } = useMsLandscapePlotQuery(params, {
  //   skip: !params,
  // });

  const [file, setFile] = useState(null);
  const fileRef = useRef(null);

  const dispatch = useDispatch();
  const mergeExposureState = (state) =>
    dispatch(actions.mergeExposure({ ...state }));
  let variableData = useSelector((state) => state.exposure.variableData) || [];

  const variableFileName =
    useSelector((state) => state.exposure.variableFileName) || [];

  function handleChange(e) {
    setFile(e.target.files[0].name);
  }

  async function Recalculate() {
    if (!fileRef) return;

    if (!fileRef.current?.files?.length) {
      mergeExposureState({
        variableFileName: '',
        variableData: null,
      });
    } else {
      const variableFileName = fileRef.current.files[0].name;
      const files = fileRef.current.files;
      const fileData = await readFile(files[0]);
      const variableData = parseMatrix(fileData);
      mergeExposureState({
        variableFileName: variableFileName,
        variableData: variableData,
      });
    }
  }

  function removeFile() {
    if (!fileRef.current?.files?.length) return;
    fileRef.current.files.value = '';
    setFile();
    fileRef.current.value = null;
    //mergeExposureState({ variableData: fileRef.current.value });
  }

  const [calculationQuery, setCalculationQuery] = useState('');
  let { data, error, isFetching } = useMsLandscapePlotQuery(calculationQuery, {
    skip: !calculationQuery,
  });
  if (data) {
    data = MsLandscape(
      data.output.cosineData,
      data.output.exposureData,
      variableData
    );
  }
  // useEffect(() => {
  //   const { study, strategy, signatureSetName } = publicForm;
  //   if (study) {
  //     setParams({
  //       study: study.value,
  //       strategy: strategy.value,
  //       signatureSetName: signatureSetName.value,
  //     });
  //   }
  // }, [publicForm]);
  useEffect(() => {
    const { study, strategy, signatureSetName, cancer, id } = state;
    if (study) {
      setCalculationQuery({
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        cancer: cancer.value,
      });
      // } else if (id) {
      //   setCalculationQuery({ userId: id });
    }
  }, [state]);

  return (
    <>
      <Form className="p-3">
        <Row className="">
          <Col lg="auto">
            <Group controlId="landscape">
              <Label>
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
              </Label>
              <div className="d-flex">
                <Form.File
                  ref={fileRef}
                  id="variableData"
                  label={file || 'Upload here (optional)'}
                  title={file || 'Upload here (optional)'}
                  onChange={(e) => {
                    handleChange(e);
                  }}
                  custom
                />
                {fileRef.current?.files?.length > 0 && (
                  <Button
                    className="ml-1"
                    size="sm"
                    title="Remove"
                    variant="danger"
                    onClick={removeFile}
                  >
                    <FontAwesomeIcon icon={faFolderMinus} size="lg" />
                  </Button>
                )}
              </div>
            </Group>
          </Col>
          <Col lg="auto" className="d-flex">
            <Button
              // disabled={source == 'user' && !id}
              className="mt-auto mb-3"
              variant="primary"
              onClick={Recalculate}
            >
              Recalculate
            </Button>
          </Col>
        </Row>
      </Form>
      <hr />
      <LoadingOverlay active={isFetching} />
      {data && !error ? (
        <Plotly
          className="w-100"
          data={data.traces}
          layout={data.layout}
          config={data.config}
        />
      ) : (
        <div className="text-center my-4">No data available</div>
      )}
    </>
  );
}
