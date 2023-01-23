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
import { useForm, Controller } from 'react-hook-form';
import { useRef } from 'react';

import './plot.scss';
import MsLandscape from '../../../controls/plotly/msLandscape/msLandscape';
import { readFile, parseMatrix } from '../../../controls/utils/utils';
import { actions } from '../../../../services/store/exposure';

const { Label, Group } = Form;
export default function MsLandscapePlot({ calculateLandscape }) {
  const publicForm = useSelector((state) => state.exposure.publicForm);
  const exposure = useSelector((state) => state.exposure);
  const { variableFile, plotPath, debugR, err, loading } = exposure.msLandscape;
  // const [params, setParams] = useState('');
  // const { data, error, isFetching } = useMsLandscapePlotQuery(params, {
  //   skip: !params,
  // });
  const fileRef = useRef(null);

  const dispatch = useDispatch();
  const mergeExposureState = (state) =>
    dispatch(actions.mergeExposure({ ...state }));
  const variableData =
    useSelector((state) => state.exposure.variableData) || [];

  const variableFileName =
    useSelector((state) => state.exposure.variableFileName) || [];

  function handleVariableData(event) {
    const variableFileName = event.target.files[0].name;
    mergeExposureState({ variableFileName });
  }

  async function Recalculate() {
    if (!fileRef) return;
    // const text = await readFile(event.target.files[0]);
    console.log(fileRef);
    console.log(fileRef.current);

    if (!fileRef.current?.files?.length) return;
    const files = fileRef.current.files;
    console.log(files);
    const fileData = await readFile(files[0]);
    console.log(fileData);
    const variableData = parseMatrix(fileData);
    console.log(variableData);
    mergeExposureState({ variableData });
  }

  function removeFile() {
    debugger;
    if (!fileRef.current?.files?.length) return;
    fileRef.current.files.value = '';
    variableFileName = '';
    variableData = [];
    mergeExposureState({ variableFileName, variableData });
  }

  console.log(variableFileName);
  console.log(variableData);
  const [calculationQuery, setCalculationQuery] = useState('');
  let { data, error, isFetching } = useMsLandscapePlotQuery(calculationQuery, {
    skip: !calculationQuery,
  });
  if (data) {
    console.log(data);
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
    const { study, strategy, signatureSetName, cancer } = publicForm;
    if (study) {
      setCalculationQuery({
        study: study.value,
        strategy: strategy.value,
        signatureSetName: signatureSetName.value,
        cancer: cancer.value,
      });
    }
  }, [publicForm]);

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
                  label={variableFileName || 'Upload here (optional)'}
                  title={variableFileName || 'Upload here (optional)'}
                  type="input"
                  onChange={(e) => handleVariableData(e)}
                  custom
                />
                {fileRef.current?.files?.length > 0 && (
                  <Button
                    className="ml-1"
                    size="sm"
                    title="Remove"
                    variant="danger"
                    disabled={loading}
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
              // disabled={source == 'user' && !projectID}
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
