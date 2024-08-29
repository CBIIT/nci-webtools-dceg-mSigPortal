import { visualizationApiSlice } from '../../../services/store/rootApi';
import { groupBy, startCase } from 'lodash';

export const publicationsSlice = visualizationApiSlice.injectEndpoints({
  endpoints: (builder) => ({
    publications: builder.query({
      query: () => ({ url: `publications` }),
      transformResponse: (data) => {
        function customCell(label = '', url = '') {
          if (url && String(url).includes('http')) {
            return (
              <a href={url} target="_blank" rel="noreferrer">
                {label}
              </a>
            );
          } else return label;
        }

        const reducer = (acc, column) => [
          ...acc,
          {
            Header:
              column === column.toLowerCase()
                ? column.length < 4
                  ? column.toUpperCase()
                  : startCase(column)
                : startCase(column),
            accessor: column,
            id: column,
            Cell: (e) => {
              // convert to hyperlinks if possible
              if (column == 'title' && e.row.values['doi']) {
                return customCell(e.value, e.row.values['doi']);
              } else if (
                column == 'softwareName' &&
                e.row.values['sourceUrl']
              ) {
                return customCell(e.value, e.row.values['sourceUrl']);
              } else {
                return customCell(e.value, e.value);
              }
            },
          },
        ];

        const columnOrder = {
          'Original Research A': [
            'diseaseOrPhenotypeOrExposure',
            'cancerType',
            'experimentalStrategy',
            'firstAuthor',
            'lastAuthor',
            'year',
            'journal',
            'bioRxivOrPubmedId',
            'title',
            'note',
            'doi',
          ],
          'Original Research B': [
            'diseaseOrPhenotypeOrExposure',
            'cancerType',
            'experimentalStrategy',
            'firstAuthor',
            'lastAuthor',
            'year',
            'bioRxivOrPubmedId',
            'title',
            'doi',
          ],
          'Review Paper': [
            'firstAuthor',
            'lastAuthor',
            'year',
            'journal',
            'bioRxivOrPubmedId',
            'title',
            'note',
            'doi',
          ],
          'Computational Methods': [
            'softwareName',
            'computationalMethod',
            'programmingLanguage',
            'firstAuthor',
            'lastAuthor',
            'year',
            'bioRxivOrPubmedId',
            'title',
            'doi',
            'sourceUrl',
          ],
        };

        const groupData = groupBy(data, (e) => e.category);
        const tableData = Object.entries(groupData).reduce(
          (obj, [category, data]) => {
            const filterData = data.map((e) => {
              const { category, id, ...rest } = e;
              return rest;
            });

            return {
              ...obj,
              [category]: {
                data: filterData,
                columns: columnOrder[category].reduce(reducer, []),
              },
            };
          },
          {}
        );

        return tableData;
      },
    }),
  }),
});

export const { usePublicationsQuery } = publicationsSlice;
