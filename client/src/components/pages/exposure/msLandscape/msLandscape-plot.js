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
import { useSelector } from 'react-redux';
import { useMsLandscapePlotQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

import './plot.scss';
import MsLandscape from '../../../controls/plotly/msLandscape/msLandscape';

const { Label, Group } = Form;
export default function MsLandscapePlot({
  calculateLandscape,
  handleVariable,
}) {
  const publicForm = useSelector((state) => state.exposure.publicForm);
  const exposure = useSelector((state) => state.exposure);
  const { variableFile, plotPath, debugR, err, loading } = exposure.msLandscape;
  // const [params, setParams] = useState('');
  // const { data, error, isFetching } = useMsLandscapePlotQuery(params, {
  //   skip: !params,
  // });

  console.log(variableFile);
  const [calculationQuery, setCalculationQuery] = useState('');
  console.log();
  let { data, error, isFetching } = useMsLandscapePlotQuery(calculationQuery, {
    skip: !calculationQuery,
  });
  console.log(data);
  if (data) {
    console.log(data);
    data = MsLandscape(data.output.cosineData, data.output.exposureData);
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

  function readFile(file) {
    let fileReader = new FileReader();
    return new Promise((resolve, reject) => {
      fileReader.onload = () => {
        resolve(fileReader.result);
      };
      fileReader.onerror = reject;
      fileReader.readAsText(file);
    });
  }

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
                  id="variableData"
                  label={variableFile || 'Upload here (optional)'}
                  title={variableFile || 'Upload here (optional)'}
                  value={''}
                  accept=""
                  onChange={(e) => handleVariable(e.target.files[0])}
                  custom
                />
                {variableFile && (
                  <Button
                    className="ml-1"
                    size="sm"
                    title="Remove"
                    variant="danger"
                    disabled={loading}
                    onClick={() => {
                      handleVariable(new File([], ''));
                    }}
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
              onClick={calculateLandscape}
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
