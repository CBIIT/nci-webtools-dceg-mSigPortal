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
        const regex4 = /^.{0,4}/;
        let showXlable = true;
        let xTitle = "";
        let yTitle = "";

        const groupByMutation = data.reduce((groups, e) => {
          let mutation;
          if (profile.label === "SBS84") {
            mutation = e.MutationType.match(regex)[1];
          } else if (profile.label === "DBS1") {
            mutation = e.MutationType.match(regex2)[0];
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

        let groupByFirstGroup = {},
          groupByMutationID = {},
          groupR = {},
          groupRDel = {},
          groupRIns = {},
          groupM = {},
          annotationIDTop = {},
          annotationIDBot = {},
          annotationIDBotArr = {},
          arrayID = [],
          arrayIDAnnotation = [],
          arrayIDAnnotationTop = [];

        if (profile.label === "ID1") {
          showXlable = false;
          xTitle = "";
          yTitle = "Number of Indels";
          console.log(Object.values(groupByMutation)[0]);

          groupByFirstGroup = Object.fromEntries(
            Object.entries(groupByMutation).slice(0, 4)
          );

          groupByMutationID = data.reduce((groups, e) => {
            let mutationID;
            mutationID = e.MutationType.match(
              e.MutationType.substring(
                e.MutationType.length - 3,
                e.MutationType.length - 2
              )
            );

            const signature = {
              mutationType: e.MutationType,
              contribution: e.Contribution,
            };
            groups[mutationID] = groups[mutationID]
              ? [...groups[mutationID], signature]
              : [signature];
            return groups;
          }, {});

          groupR = groupByMutationID["R"].reduce((r, a) => {
            let m;
            m = a.mutationType.match(a.mutationType.substr(2, 3));
            const s = {
              mutationType: a.mutationType,
              contribution: a.contribution,
            };
            r[m] = r[m] ? [...r[m], a] : [s];
            return r;
          }, {});

          groupRDel = groupR["Del"].reduce((r, a) => {
            let m;
            m = a.mutationType.match(a.mutationType.substr(0, 7));
            const s = {
              mutationType: a.mutationType,
              contribution: a.contribution,
            };
            r[m] = r[m] ? [...r[m], a] : [s];
            return r;
          }, {});

          groupRIns = groupR["Ins"].reduce((r, a) => {
            let m;
            m = a.mutationType.match(a.mutationType.substr(0, 7));
            const s = {
              mutationType: a.mutationType,
              contribution: a.contribution,
            };
            r[m] = r[m] ? [...r[m], a] : [s];
            return r;
          }, {});

          groupM = groupByMutationID["M"].reduce((r, a) => {
            let m;
            m = a.mutationType.match(a.mutationType.substr(0, 7));
            const s = {
              mutationType: a.mutationType,
              contribution: a.contribution,
            };
            r[m] = r[m] ? [...r[m], a] : [s];
            return r;
          }, {});
          const arrayID1 = Object.keys(groupByFirstGroup).map(function (key) {
            return groupByFirstGroup[key];
          });
          const arrayID2 = Object.keys(groupRDel).map(function (key) {
            return groupRDel[key];
          });
          const arrayID3 = Object.keys(groupRIns).map(function (key) {
            return groupRIns[key];
          });
          const arrayID4 = Object.keys(groupM).map(function (key) {
            return groupM[key];
          });

          arrayID = [...arrayID1, ...arrayID2, ...arrayID3, ...arrayID4];

          Object.values(arrayID).map((group) =>
            group.map((e) => arrayIDAnnotation.push(e.mutationType))
          );
          console.log(arrayIDAnnotation);

          console.log(groupByFirstGroup);
          console.log(arrayID1);
          console.log(arrayID2);
          console.log(arrayID3);
          console.log(arrayID4);
          console.log(arrayID);
          // annotationsIDTop = Object.entries(groupByMutationID).map(
          //   ([mutation, signatures]) => ({
          //     xref: "x",
          //     yref: "paper",
          //     x: signatures.map((e) => e.mutationType)[signatures.length / 2],
          //     xanchor: "bottom",
          //     y: 1.05,
          //     yanchor: "bottom",
          //     text: "<b>" + mutation + "</b>",
          //     showarrow: false,
          //     font: {
          //       size: 14,
          //     },
          //     align: "center",
          //   })
          // );
          annotationIDTop = Object.entries(arrayID).map(
            ([mutation, index]) => ({
              x: index.map((e) =>
                index < 4 ? e.mutationType[index.length / 2] : e.mutationType[0]
              ),
              y: mutation,
            })
          );

          console.log(annotationIDTop);

          annotationIDBot = arrayIDAnnotation.map((num) => {
            const result = {
              xref: "x",
              yref: "paper",
              x: num,
              xanchor: "bottom",
              y: -0.08,
              yanchor: "bottom",
              text: "<b>" + num.substr(num.length - 1, num.length) + "</b>",
              showarrow: false,
              font: {
                size: 14,
              },
              align: "center",
            };
            return result;
          });
          console.log(annotationIDBot);
        } else if (profile.label === "SBS84") {
          showXlable = true;
          xTitle = "Substitution";
          yTitle = "Mutation Probability";
        } else {
          showXlable = true;
          xTitle = "Double substitution";
          yTitle = "Mutation Probability";
        }

        console.log(annotationIDBot);

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
          "1:Del:C": "#FBBD6F",
          "1:Del:T": "#FE8002",
          "1:Ins:C": "#AEDD8A",
          "1:Ins:T": "#35A12E",
          "2:Del:R": "#FCC9B4",
          "3:Del:R": "#FB8969",
          "4:Del:R": "#F04432",
          "5:Del:R": "#BB1A1A",
          "2:Ins:R": "#CFDFF0",
          "3:Ins:R": "#93C3DE",
          "4:Ins:R": "#4B97C7",
          "5:Ins:R": "#1863AA",
          "2:Del:M": "#E1E1EE",
          "3:Del:M": "#B5B5D6",
          "4:Del:M": "#8482BC",
          "5:Del:M": "#62409A",
        };

        const annotationColors = {
          "1:Del:C": "black",
          "1:Del:T": "white",
          "1:Ins:C": "black",
          "1:Ins:T": "white",
          "2:Del:R": "black",
          "3:Del:R": "black",
          "4:Del:R": "black",
          "5:Del:R": "white",
          "2:Ins:R": "black",
          "3:Ins:R": "black",
          "4:Ins:R": "black",
          "5:Ins:R": "white",
          "2:Del:M": "blacl",
          "3:Del:M": "black",
          "4:Del:M": "black",
          "5:Del:M": "white",
        };

        const annText = (mutation, res) => {
          if (profile.label === "SBS84") {
            res = "<b>" + mutation + "</b>";
          } else if (profile.label === "DBS1") {
            res = "<b>" + mutation + "NN</b>";
          } else if (profile.label === "ID1") {
            res = "<b>" + mutation.substring(mutation.length - 1) + "</b>";
          }
          return res;
        };

        const traces1 = Object.entries(groupByMutation).map(
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

        const traces2 = Object.entries(groupByFirstGroup).map(
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
        const traces3 = Object.entries(groupRDel).map(
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

        const traces4 = Object.entries(groupRIns).map(
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

        const traces5 = Object.entries(groupM).map(
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

        let traces = [];

        const shapes1 = Object.entries(groupByMutation).map(
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

        const shapes2 = Object.entries(groupByMutation).map(
          ([mutation, signatures]) => ({
            type: "rect",
            xref: "x",
            yref: "paper",
            x0: signatures.map((e) => e.mutationType)[0],
            y0: -0.01,
            x1: signatures.map((e) => e.mutationType)[signatures.length - 1],
            y1: -0.04,
            fillcolor: colors[mutation],
            line: {
              width: 0,
            },
            //layer: "below",
          })
        );

        const annotationsOthers = Object.entries(groupByMutation).map(
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
              //color: colors[mutation],
              size: 18,
            },
            align: "center",
          })
        );

        console.log(groupByMutation);

        const annotationsID = Object.entries(groupByMutation).map(
          ([mutation, signatures]) => ({
            xref: "x",
            yref: "paper",
            x: signatures.map((e) => e.mutationType)[signatures.length / 2],
            xanchor: "bottom",
            y: 1,
            yanchor: "bottom",
            text: annText(mutation),
            showarrow: false,
            font: {
              color: annotationColors[mutation],
              size: 14,
            },
            align: "center",
          })
        );

        console.log(annotationsID);
        let annotations = [];
        let shapes = [...shapes1];

        if (profile.label === "ID1") {
          shapes = [...shapes1, ...shapes2];
          annotations = [...annotationsID, ...annotationIDBot];
          traces = [...traces2, ...traces3, ...traces4, ...traces5];
        } else {
          annotations = [...annotationsOthers];
          traces = [...traces1];
        }

        console.log(traces);
        const layout = {
          xaxis: {
            title: xTitle,
            showline: true,
            showticklabels: showXlable,
            tickangle: -90,
            tickfont: {
              size: 10,
            },
          },
          yaxis: {
            title: yTitle,
            //tickformat: ".1%",
            autorange: true,
          },

          shapes: shapes,
          annotations: annotations,
        };

        return { data: traces, layout, config: {} };
      } else return { data: [], layout: {}, config: {} };
    } catch (error) {
      return { data: [], layout: {}, config: {} };
    }
  },
});
