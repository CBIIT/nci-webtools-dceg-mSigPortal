import { useEffect, useState, useMemo } from 'react';
import { useForm } from 'react-hook-form';
import Select from '../../../controls/select/selectForm';
import { Container, Form, Row, Col } from 'react-bootstrap';
import Plotly from '../../../controls/plotly/plot/plot';
import { useSelector, useDispatch } from 'react-redux';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Description from '../../../controls/description/description';
import { useProfilerSummaryQuery } from './apiSlice';
import { actions } from '../../../../services/store/visualization';

export default function ProfilerSummary() {
  const store = useSelector((state) => state.visualization);
  const publicForm = store.publicForm;
  const { matrixData, source, id } = store.main;
  const { filter } = store.profilerSummary;

  const dispatch = useDispatch();
  const mergeState = (state) =>
    dispatch(actions.mergeVisualization({ profilerSummary: state }));
  const [params, setParams] = useState();

  const { data, error, isFetching } = useProfilerSummaryQuery(params, {
    skip: !params,
  });

  const { control } = useForm();

  const filterOptions = useMemo(() => {
    const filters = matrixData.length
      ? [...new Set(matrixData.map((d) => d.filter))].map((e) => ({
          label: e ? e : 'N/A',
          value: e,
        }))
      : [];
    return filters.length > 1
      ? [...filters, { label: 'All', value: false }]
      : filters;
  }, [matrixData]);

  // automatically select first filter option
  useEffect(() => {
    if (filterOptions.length) mergeState({ filter: filterOptions[0] });
  }, [filterOptions]);

  // query new summary plot when data is recieved or filter is changed
  useEffect(() => {
    if (matrixData.length) {
      const params = {
        study: publicForm.study.value,
        cancer: publicForm.cancer.value,
        strategy: publicForm.strategy.value,
        ...(source == 'user' && { userId: id }),
        ...(filter?.value !== false && { filter: filter.value }),
      };
      if (source == 'public' || (source == 'user' && filter?.label)) {
        setParams(params);
      }
    }
  }, [matrixData, filter]);

  return (
    <div className="bg-white border rounded" style={{ minHeight: '500px' }}>
      <div className="p-3">
        <b>Number of Mutations Per Sample with Regard to Mutational Profile</b>
        <Description
          className="m-0"
          less="This plot illustrates the number of mutations in each tumor sample from [Cancer Type] in the selected [Study]."
          more="On the y-axis is the number of mutations in log base 10 scale, and on the x-axis is the sample index for each sample of the selected cancer type (sorted by number of mutations in ascending order). The different colored lines represent different mutational profiles (SBS= single-base substitution, DBS= doublet-base substitution, ID=indel)."
        />
      </div>
      <hr />
      {source == 'user' && filterOptions.length && (
        <>
          <Container fluid>
            <Form>
              <Row>
                <Col sm="auto">
                  <Select
                    className="mb-2"
                    name="filter"
                    label="Filter"
                    value={filter || filterOptions[0]}
                    disabled={isFetching}
                    options={filterOptions}
                    control={control}
                    onChange={(filter) => mergeState({ filter })}
                  />
                </Col>
              </Row>
            </Form>
          </Container>
          <hr />
        </>
      )}
      <LoadingOverlay active={isFetching} />
      {data && (
        <Plotly
          data={data.traces}
          layout={data.layout}
          config={data.config}
          filename={
            source == 'public'
              ? `${publicForm.study.value}_profiler-summary`
              : 'profiler-summary'
          }
        />
      )}
      {error && (
        <p className="text-center">
          An error has occured. Please check your inputs and try again.
        </p>
      )}
    </div>
  );
}
