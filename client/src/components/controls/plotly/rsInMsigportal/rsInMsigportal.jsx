import { groupBy } from 'lodash';

export default function RsInMsigportal(rawData) {
  const COL_START = 0.02; // left edge of the first column
  const COL_PITCH = 0.138; // horizontal distance between adjacent column left edges
  const PIE_W = 0.096; // width of each pie's x-domain
  const ROW_H = 0.25; // vertical distance between rows
  const PIE_PAD = 0.02; // vertical padding inside each row
  const TITLE_PAD = 0.015; // gap between a pie's top and its profile title
  const ROWS_PER_COL = 4; // rows before wrapping into the next column

  const colX = (g) => [
    COL_START + g * COL_PITCH,
    COL_START + g * COL_PITCH + PIE_W,
  ];
  const colCenter = (g) => COL_START + g * COL_PITCH + PIE_W / 2;
  const rowY = (r) => [1 - (r + 1) * ROW_H + PIE_PAD, 1 - r * ROW_H - PIE_PAD];
  const titleY = (r) => 1 - r * ROW_H - TITLE_PAD;

  // ---- species layout config (display order, left -> right) ----
  // `match` is the binomial name so it is robust to genome-assembly suffix
  const SPECIES = [
    {
      match: 'Homo sapiens',
      profileOrder: [
        'SBS96',
        'SBS192',
        'SBS288',
        'SBS1536',
        'DBS78',
        'ID29',
        'ID83',
        'CN48',
        'RS32',
        'RNA192',
      ],
      startCol: 0,
      cols: 3,
      titlePad: 8,
    },
    {
      match: 'Mus musculus',
      profileOrder: ['SBS96', 'DBS78', 'ID83'],
      startCol: 3,
      cols: 1,
      titlePad: 7,
    },
    {
      match: 'Rattus',
      profileOrder: ['SBS96', 'DBS78'],
      startCol: 4,
      cols: 1,
      titlePad: 7,
    },
    {
      match: 'Gallus',
      profileOrder: ['SBS96'],
      startCol: 5,
      cols: 1,
      titlePad: 7,
    },
    {
      match: 'Caenorhabditis',
      profileOrder: ['SBS96'],
      startCol: 6,
      cols: 1,
      titlePad: 7,
    },
  ];

  const HOVER =
    '<b>%{label}</b> <br>%{percent} </br> %{value}  <extra></extra>';

  const bySpecies = groupBy(rawData, (d) => `${d.species}`);
  const setNames = [...new Set(rawData.map((d) => d.signatureSetName))];
  const colors = {};
  const usedColors = new Set();
  setNames.forEach((name) => {
    let hex;
    do {
      hex =
        '#' +
        ('000000' + Math.floor(Math.random() * 0x1000000).toString(16)).slice(
          -6
        );
    } while (usedColors.has(hex));
    usedColors.add(hex);
    colors[name] = hex;
  });

  const traces = [];
  const annotations = [];
  const shapes = [];

  SPECIES.forEach((cfg, si) => {
    // the actual species key may carry an assembly suffix; match on the binomial
    const speciesKey = Object.keys(bySpecies).find((k) =>
      k.includes(cfg.match)
    );
    const endCol = cfg.startCol + cfg.cols - 1;

    if (speciesKey) {
      const byProfile = groupBy(
        bySpecies[speciesKey],
        (d) => `${d.profile}${d.matrix}`
      );
      const order =
        cfg.profileOrder && cfg.profileOrder.length
          ? cfg.profileOrder
          : Object.keys(byProfile);

      order.forEach((profileKey, idx) => {
        const rows = byProfile[profileKey];
        if (!rows || !rows.length) return;

        const g = cfg.startCol + Math.floor(idx / ROWS_PER_COL);
        const r = idx % ROWS_PER_COL;

        traces.push({
          type: 'pie',
          marker: {
            color: rows.map((e) => colors[e.signatureSetName]),
            line: { color: 'black', width: 1 },
          },
          textposition: 'inside',
          labels: rows.map((e) => e.signatureSetName),
          values: rows.map((e) => parseInt(e.count, 10)),
          texttemplate: '%{value}',
          direction: 'clockwise',
          name: profileKey,
          domain: { x: colX(g), y: rowY(r) },
          hovertemplate: HOVER,
        });

        annotations.push({
          xref: 'paper',
          yref: 'paper',
          xanchor: 'center',
          yanchor: 'bottom',
          showarrow: false,
          text: profileKey.padStart(cfg.titlePad, ' '),
          align: 'center',
          font: { weight: 'bold' },
          x: colCenter(g),
          y: titleY(r),
        });
      });

      // species header centered over its column span (assembly on its own line)
      annotations.push({
        xref: 'paper',
        yref: 'paper',
        xanchor: 'center',
        yanchor: 'bottom',
        showarrow: false,
        text: speciesKey.replace(/\s*\(([^)]+)\)/, '<br>($1)'),
        font: { size: 14, weight: 'bold' },
        x: (colCenter(cfg.startCol) + colCenter(endCol)) / 2,
        y: 1.02,
      });
    }

    // vertical separator to the right of this species (not after the last one)
    if (si < SPECIES.length - 1) {
      const sepX = (colX(endCol)[1] + colX(endCol + 1)[0]) / 2;
      shapes.push({
        type: 'line',
        xref: 'paper',
        yref: 'paper',
        x0: sepX,
        y0: 0,
        x1: sepX,
        y1: 1,
        line: { color: 'lightgray', width: 2, dash: 'solid' },
      });
    }
  });

  const layout = {
    hoverlabel: { bgcolor: '#FFF' },
    height: 920,
    width: 1600,
    autosize: false,
    margin: {
      l: 0,
      r: 50,
      t: 100,
      b: 50,
    },
    legend: {
      title: {
        text: '\t <b>Signature Set Name</b>',
        font: {
          family: 'Times New Roman',
          size: 17,
        },
      },
      x: 1,
      xanchor: 'right',
      y: 0,
    },
    annotations,
    shapes,
  };
  const config = {
    responsive: true,
    displayModeBar: true,
  };
  return { traces, layout, config };
}
