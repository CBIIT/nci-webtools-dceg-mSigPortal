import { useEffect, useState } from 'react';
import { Button, Form, Row, Col } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { useForm, Controller } from 'react-hook-form';

import Plotly from '../../../controls/plotly/plot/plot';
import { actions } from '../../../../services/store/visualization';
import { useMpeaScatterQuery, useMpeaBarQuery } from './apiSlice';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

export default function MutPatternPlot() {
  const store = useSelector((state) => state.visualization);
  const { study, cancer, strategy } = store.publicForm;
  const { source, id, matrixList } = store.main;
  const { proportion, pattern } = store.mutationalPattern;

  const [scatterParams, setScatterParams] = useState('');

  const dispatch = useDispatch();
  const mergeState = (state) =>
    dispatch(actions.mergeVisualization({ mutationalPattern: state }));

  const {
    control,
    handleSubmit,
    formState: { errors },
  } = useForm({
    defaultValues: {
      proportion: proportion || 0.8,
      pattern: pattern || 'NCG>NTG',
    },
  });

  const {
    data: scatterData,
    error: scatterError,
    isFetching: fetchingScatter,
  } = useMpeaScatterQuery(scatterParams, {
    skip: !scatterParams,
  });

  const [barParams, setBarParams] = useState('');
  const {
    data: patternData,
    error: patternError,
    isFetching: fetchingPattern,
  } = useMpeaBarQuery(barParams, {
    skip: !barParams,
  });

  // get data on form change
  useEffect(() => {
    if (pattern) {
      setScatterParams({
        profile: 'SBS',
        matrix: '96',
        pattern: pattern,
        ...(source == 'public' && {
          study: study.value,
          cancer: cancer.value,
          strategy: strategy.value,
        }),
        ...(source == 'user' && { userId: id }),
      });
    }
  }, [study, proportion, pattern]);

  useEffect(() => {
    if (source == 'public') {
      if (study && proportion) {
        setBarParams({
          study: study.value,
          proportion: proportion,
        });
      }
    } else {
      if (proportion && pattern) {
        setBarParams({
          pattern,
          proportion,
          matrixFile: matrixList.filter(
            (e) => e.profile == 'SBS' && e.matrix == '96'
          )[0].Path,
          userId: id,
        });
      }
    }
  }, [proportion, pattern]);

  async function onSubmit(data) {
    const { pattern, proportion } = data;
    mergeState({ pattern, proportion: parseFloat(proportion) });
  }

  return (
    <div>
      <Form className="p-3" onSubmit={handleSubmit(onSubmit)}>
        <Row>
          <Col lg="auto">
            <Form.Group
              controlId="minimum"
              title="Minimal Proportion mutations within Each Mutational Pattern"
            >
              <Form.Label>Minimal Proportion (0-1)</Form.Label>
              <Controller
                name="proportion"
                control={control}
                rules={{ required: true, min: 0, max: 1 }}
                render={({ field }) => (
                  <Form.Control
                    {...field}
                    type="number"
                    min="0"
                    max="1"
                    step="0.01"
                    placeholder="Ex. 0.8"
                    isInvalid={errors.proportion}
                    disabled={fetchingPattern || fetchingScatter}
                  />
                )}
              />
              <Form.Control.Feedback type="invalid">
                {errors.proportion?.type == 'required'
                  ? 'Minimal Proportion is required'
                  : 'Enter a value between 0 and 1'}
              </Form.Control.Feedback>
            </Form.Group>
          </Col>
          <Col lg="auto">
            <Form.Group controlId="pattern">
              <Form.Label>Mutational Pattern</Form.Label>
              <Controller
                name="pattern"
                control={control}
                rules={{ required: true }}
                render={({ field }) => (
                  <Form.Control
                    {...field}
                    placeholder="Ex. NCG>NTG"
                    isInvalid={errors.pattern}
                    disabled={fetchingPattern || fetchingScatter}
                  />
                )}
              />
              <Form.Control.Feedback type="invalid">
                {errors.pattern?.type == 'required' && 'Enter a valid pattern'}
              </Form.Control.Feedback>
            </Form.Group>
          </Col>
          <Col lg="auto" className="d-flex">
            <Button
              className="mt-auto mb-3"
              variant="primary"
              type="submit"
              disabled={fetchingPattern || fetchingScatter}
            >
              Calculate
            </Button>
          </Col>
        </Row>
      </Form>
      <hr />
      <div id="barchart">
        {fetchingPattern && (
          <div style={{ minHeight: '150px' }}>
            <LoadingOverlay
              active={fetchingPattern}
              showIndicator={true}
              content="Loading Frequency Plot"
            />
          </div>
        )}
        {patternData && patternData?.traces && !patternError && (
          <>
            <Plotly
              className="w-100"
              data={patternData.traces}
              layout={patternData.layout}
              config={patternData.config}
            />
            <p className="p-3">
              This plot illustrates the frequency by count of each mutational
              pattern in the given study and cancer type or input dataset. The
              y-axis is the frequency of each mutational pattern across all
              samples, and the x-axis includes each of the mutational patterns
              present in the study and cancer type that meet the criteria for
              the minimal proportion of mutations within each mutational
              pattern.
            </p>
          </>
        )}
        {patternData && !patternData?.traces && !fetchingPattern && (
          <div className="p-3 text-center">
            <b>Frequency of Mutational Pattern</b>
            <p>
              No mutational pattern with proportion of mutations large than{' '}
              {proportion}. Try a lower proportion value.
            </p>
          </div>
        )}
        {patternError && (
          <div className="p-3 text-center">
            <b>Frequency of Mutational Pattern</b>
            <p>An error has occured. Please verify your input.</p>
          </div>
        )}
      </div>
      <hr />
      <div id="scatter">
        {fetchingScatter && (
          <div style={{ minHeight: '150px' }}>
            <LoadingOverlay active={fetchingScatter} />
          </div>
        )}
        {scatterData && !scatterError && (
          <>
            <Plotly
              className="w-100"
              data={scatterData.traces}
              layout={scatterData.layout}
              config={scatterData.config}
            />
            <p className="p-3">
              This plot illustrates the mutational pattern context entered
              compared to other contexts with the same SBS mutation for each
              sample. On the y-axis is the other contexts, and on the x-axis is
              the specific mutational pattern context input for the enrichment
              analysis. For some studies including multiple cancer types (such
              as TCGA PanCancer), different colors will be used.
            </p>
          </>
        )}
        {patternError && (
          <div className="p-3 text-center">
            <b>Proportion of Mutational Pattern</b>
            <p>An error has occured. Please verify your input.</p>
          </div>
        )}
      </div>
    </div>
  );
}
