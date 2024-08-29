import { useEffect, useState, useMemo } from 'react';
import { useForm } from 'react-hook-form';
import Select from '../../../controls/select/selectHookForm';
import { Container, Form, Row, Col } from 'react-bootstrap';
import Plotly from '../../../controls/plotly/plot/plot';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';
import Description from '../../../controls/description/description';
import { useProfilerSummaryQuery } from './apiSlice';
import { useSeqmatrixOptionsQuery } from '../../../../services/store/rootApi';

export default function ProfilerSummary({ state }) {
  const { source, study, cancer, strategy, id } = state;
  const [params, setParams] = useState();

  const { data: options } = useSeqmatrixOptionsQuery(
    {
      ...(source == 'public'
        ? { study: study.value, cancer: cancer.value, strategy: strategy.value }
        : { userId: id }),
    },
    { skip: source == 'user' ? !id : !study }
  );
  const { data, error, isFetching } = useProfilerSummaryQuery(params, {
    skip: !params,
  });

  const { control, watch, setValue } = useForm({
    defaultValues: { filter: '' },
  });
  const { filter } = watch();

  const filterOptions = useMemo(() => {
    const filters = options
      ? [...new Set(options.map((d) => d.filter))].map((e) => ({
          label: e ? e : 'N/A',
          value: e,
        }))
      : [];
    return filters.length > 1
      ? [...filters, { label: 'All', value: false }]
      : filters;
  }, [options]);

  // automatically select first filter option
  useEffect(() => {
    if (filterOptions.length && !filter) setValue('filter', filterOptions[0]);
  }, [filterOptions]);

  // query new summary plot when data is received or filter is changed
  useEffect(() => {
    if (options) {
      const params = {
        study: study.value,
        cancer: cancer.value,
        strategy: strategy.value,
        ...(source == 'user' && { userId: id }),
        ...(filter?.value !== false && { filter: filter.value }),
      };

      setParams(params);
    }
  }, [options, filter]);

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
              ? `${state.study.value}_profiler-summary`
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
