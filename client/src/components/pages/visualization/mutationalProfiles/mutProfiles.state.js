import { atom, selector } from 'recoil';
import axios from 'axios';
import { visualizationState } from '../visualization.state';

export const defaultFormState = {
  profile: {},
  // filtered: [],
  // selectName: '',
  // selectProfile: '',
  // selectMatrix: '',
  // selectFilter: '',
};

export const formState = atom({
  key: 'visualization.mutProfiles.formState',
  default: defaultFormState,
});

// const defaultOptions = {
//   nameOptions: [],
//   profileOptions: [],
//   matrixOptions: [],
//   filterOptions: [],
// };

// export const formOptions = selector({
//   key: 'visualization.mutProfiles.formOptions',
//   get: async ({ get }) => {
//     const { svgList } = get(visualizationState);
//     try {
//       const { data } = await axios.post('api/visualizationWrapper', {
//         fn: 'getTreeLeaf',
//       });

//       return data.output;
//       return [];
//     } catch (error) {
//       return null;
//     }
//   },
// });

export const getPlot = selector({
  key: 'visualization.mutProfiles.plot',
  get: async ({ get }) => {
    try {
      const { profile } = get(formState);

      if (Object.keys(profile).length) {
        const { data } = await axios.post('api/querySignature', {
          filter: { ...profile.value }, // filter given an object of key:values you want
          properties: ['MutationType', 'Contribution'], // return objects containing these properties
        });

        const regex = /\[(.*)\]/;
        const groupByMutation = data.reduce((groups, e) => {
          const mutation = e.MutationType.match(regex)[1];
          const signature = {
            mutationType: e.MutationType,
            contribution: e.Contribution,
          };
          groups[mutation] = groups[mutation]
            ? [...groups[mutation], signature]
            : [signature];
          return groups;
        }, {});

        const colors = {
          'C>A': 'blue',
          'C>G': 'black',
          'C>T': 'red',
          'T>A': 'grey',
          'T>C': 'green',
          'T>G': 'pink',
        };

        const traces = Object.entries(groupByMutation).map(
          ([mutation, signatures]) => ({
            name: mutation,
            type: 'bar',
            marker: { color: colors[mutation] },
            x: signatures.map((e) => e.mutationType),
            y: signatures.map((e) => e.contribution),
            hoverinfo: 'x+y',
            showlegend: false,
          })
        );

        const shapes = Object.entries(groupByMutation).map(
          ([mutation, signatures]) => ({
            type: 'rect',
            xref: 'x',
            yref: 'paper',
            x0: signatures.map((e) => e.mutationType)[0],
            y0: 1.03,
            x1: signatures.map((e) => e.mutationType)[signatures.length - 1],
            y1: 1,
            fillcolor: colors[mutation],
            line: {
              width: 0,
            },
          })
        );

        const annotations = Object.entries(groupByMutation).map(
          ([mutation, signatures]) => ({
            xref: 'x',
            yref: 'paper',
            x: signatures.map((e) => e.mutationType)[
              Math.round(signatures.length / 2)
            ],
            xanchor: 'bottom',
            y: 1.05,
            yanchor: 'bottom',
            text: '<b>' + mutation + '</b>',
            showarrow: false,
            font: {
              // color: colors[mutation],
              size: 18,
            },
            align: 'center',
          })
        );

        const layout = {
          xaxis: {
            title: 'Substitution',
            showline: true,
            showticklabels: true,
            tickangle: -90,
            tickfont: {
              size: 10,
            },
          },
          yaxis: {
            title: 'Mutation probability',
            //tickformat: ".1%",
            autorange: true,
          },

          shapes: shapes,
          annotations: annotations,
        };

        return { data: [...traces], layout, config: {} };
      } else return { data: [], layout: {}, config: {} };
    } catch (error) {
      return { data: [], layout: {}, config: {} };
    }
  },
});
