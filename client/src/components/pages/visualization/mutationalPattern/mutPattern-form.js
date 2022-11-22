import { Button, Form, Row, Col } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { useForm, Controller } from 'react-hook-form';
import { NavHashLink } from 'react-router-hash-link';
import Description from '../../../controls/description/description';
import { actions as visualizationActions } from '../../../../services/store/visualization';

const actions = { ...visualizationActions };

export default function MutationalPatternForm() {
  const store = useSelector((state) => state.visualization);
  const { proportion, pattern } = store.mutationalPattern;

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

  async function onSubmit(data) {
    const { pattern, proportion } = data;
    mergeState({ pattern, proportion: parseFloat(proportion) });
  }

  return (
    <div>
      <div className="p-3">
        <Description
          less={
            <p>
              The aim of the mutational pattern enrichment analysis is to
              determine frequency and enrichment of different types of
              mutational patterns. For more information about mutational pattern
              enrichment, click <NavHashLink to="/faq#mpea">here</NavHashLink>.
            </p>
          }
          more={
            <>
              <p>
                <i>Minimal Proportion:</i> For the “Frequency of Mutational
                Pattern” plot, set the minimal proportion of mutational patterns
                identified in all samples from selected or input study. A
                slightly high proportion, such as 0.5, is suggested.
              </p>
              <p>
                <i>Mutational Pattern:</i> For the enrichment plot of
                “Proportion of Mutational Pattern Context Compared to Other
                Contexts with the same SBS Mutation”, select the mutational
                pattern to identify the enrichment of specific mutation context
                as suggested from “Frequency of Mutational Pattern”. The
                mutational pattern supports common nucleotide symbols. Click
                here for common nucleotide symbol information.
              </p>
            </>
          }
        />
      </div>
      <hr />
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
                  />
                )}
              />

              <Form.Control.Feedback type="invalid">
                {errors.pattern?.type == 'required' && 'Enter a valid pattern'}
              </Form.Control.Feedback>
            </Form.Group>
          </Col>
          <Col lg="auto" className="d-flex">
            <Button className="mt-auto mb-3" variant="primary" type="submit">
              Calculate
            </Button>
          </Col>
        </Row>
      </Form>
    </div>
  );
}
