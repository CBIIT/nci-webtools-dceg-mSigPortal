import { atom, selector } from "recoil";
import axios from "axios";
import { visualizationState } from "../visualization.state";

export const defaultFormState = {
  profile: {},
  // filtered: [],
  // selectName: '',
  // selectProfile: '',
  // selectMatrix: '',
  // selectFilter: '',
};

export const formState = atom({
  key: "visualization.mutProfiles.formState",
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
  key: "visualization.mutProfiles.plot",
  get: async ({ get }) => {
    try {
      const { profile } = get(formState);

      if (Object.keys(profile).length) {
        const { data } = await axios.post("api/querySignature", {
          filter: { ...profile.value }, // filter given an object of key:values you want
          properties: ["MutationType", "Contribution"], // return objects containing these properties
        });

        console.log(profile);

        console.log(data);
        const regex = /\[(.*)\]/;
        const regex2 = /^.{0,3}/;
        const regex3 = /^.{0,7}/;

        const groupByMutation = data.reduce((groups, e) => {
          let mutation;
          if (profile.label === "SBS84") {
            mutation = e.MutationType.match(regex)[1];
          } else if (profile.label === "DBS1") {
            mutation = e.MutationType.match(e.MutationType.substring(0, 3));
          } else if (profile.label === "ID1") {
            mutation = e.MutationType.match(regex3)[0];
          }

          const signature = {
            mutationType: e.MutationType,
            contribution: e.Contribution,
          };
          groups[mutation] = groups[mutation]
            ? [...groups[mutation], signature]
            : [signature];
          return groups;
        }, {});

        console.log(groupByMutation);

        const colors = {
          "C>A": "#03BCEE",
          "C>G": "black",
          "C>T": "#E32926",
          "T>A": "#CAC9C9",
          "T>C": "#A1CE63",
          "T>G": "#EBC6C4",
          "AC>": "#09BCED",
          "AT>": "#0266CA",
          "CC>": "#9FCE62",
          "CG>": "#006501",
          "CT>": "#FF9898",
          "GC>": "#E22925",
          "TA>": "#FEB065",
          "TC>": "#FD8000",
          "TG>": "#CB98FD",
          "TT>": "#4C0299",
        };

        const annText = (mutation, res) => {
          if (profile.label === "SBS84") {
            res = "<b>" + mutation + "</b>";
          } else if (profile.label === "DBS1") {
            res = "<b>" + mutation + "NN</b>";
          } else if (profile.label === "ID1") {
            res = "<b>" + mutation + "</b>";
          }
          return res;
        };

        const traces = Object.entries(groupByMutation).map(
          ([mutation, signatures]) => ({
            name: mutation,
            type: "bar",
            marker: { color: colors[mutation] },
            x: signatures.map((e) => e.mutationType),
            y: signatures.map((e) => e.contribution),
            hoverinfo: "x+y",
            showlegend: false,
          })
        );

        const shapes = Object.entries(groupByMutation).map(
          ([mutation, signatures]) => ({
            type: "rect",
            xref: "x",
            yref: "paper",
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
            xref: "x",
            yref: "paper",
            x: signatures.map((e) => e.mutationType)[
              Math.round(signatures.length / 2)
            ],
            xanchor: "bottom",
            y: 1.05,
            yanchor: "bottom",
            text: annText(mutation),
            showarrow: false,
            font: {
              color: colors[mutation],
              size: 18,
            },
            align: "center",
          })
        );

        const layout = {
          xaxis: {
            title: "Substitution",
            showline: true,
            showticklabels: true,
            tickangle: -90,
            tickfont: {
              size: 10,
            },
          },
          yaxis: {
            title: "Mutation probability",
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
