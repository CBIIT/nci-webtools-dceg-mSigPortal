import React, { useEffect } from 'react';
import { HashRouter as Router, Route } from 'react-router-dom';
import { Navbar } from './controls/navbar/navbar';
import Home from './pages/home/home';
import About from './pages/about/about';
import Visualization from './pages/visualization/visualization';
import Catalog from './pages/catalog/catalog';
import Exposure from './pages/exposure/exposure';
import Refitting from './pages/refitting/refitting';
import Association from './pages/association/association';
import Publications from './pages/publications/publications';
import Faq from './pages/faq/faq';
import { ErrorModal } from './controls/error-modal/error-modal';
import { SuccessModal } from './controls/success-modal/success-modal';
import { useSelector, useDispatch } from 'react-redux';
import { actions } from '../services/store/publications';

export default function App() {
  const links = [
    {
      route: '/catalog',
      action: 'Catalog',
      title: 'Catalog',
      cardTitle: 'Signature Catalog',
      cardText: 'Signature Catalog',
      description:
        'Systematically explore any reference or update to date published signatures with different profiles, version and etiology (endogenous vs Exogenous). Intergratively explore the landscape of signature exposure in different genomic studies, including TCGA, PCAWG, and our Sherlock-Lung study.',
      image: 'assets/images/explore.png',
      navIndex: 0,
      color: '#2c71dd', // blue
    },

    {
      route: '/visualization',
      action: 'Visualization',
      title: 'Visualization',
      cardTitle: 'Signature Visualization',
      cardText: 'Visualize mutational profiles',
      description:
        'Interactively and comprehensively visualize mutation signature in both sample and study level, including different type and level of mutational profiles (SBS/INDEL/DBS/SV/CNV), PCA components and different mutational feature (kataegis mutation, mutation quality, drive gene mutation etc).',
      image: 'assets/images/visualize.png',
      navIndex: 1,
      color: '#fc8701', // orange
    },
    {
      route: '/exploration',
      action: 'Exploration',
      title: 'Exploration',
      cardTitle: 'Signature Exploration',
      cardText: 'Signature Exploration',
      description:
        'Systematically explore any reference or update to date published signatures with different profiles, version and etiology (endogenous vs Exogenous). Intergratively explore the landscape of signature exposure in different genomic studies, including TCGA, PCAWG, and our Sherlock-Lung study.',
      image: 'assets/images/explore.png',
      navIndex: 2,
      color: '#2c71dd', // blue
    },
    // {
    //   route: '/refitting',
    //   action: 'Refitting',
    //   title: 'Refitting',
    //   cardTitle: 'Signature Refitting',
    //   cardText: 'Signature Refitting',
    //   description:
    //     'Comprehensively evaluate the accuracy of mutational signature based on different statistical variables (Cosine similarity, BIC, L2 norms etc) and re-decompsite signatures using different algorithms (SigProfiler, deconstructsig, bootstrapping method).',
    //   image: 'assets/images/refitting.png',
    //   navIndex: 2,
    //   color: '#689f39', // green
    // },
    {
      route: '/association',
      action: 'Association',
      title: 'Association',
      cardTitle: 'Signature Association',
      cardText: 'Signature Association',
      description:
        'Systematically analyze and visualize the association between mutational signature exposure and genomic or epigenomic features or other sample based variables (such as clinical data) in different genomic studies.',
      image: 'assets/images/association.png',
      navIndex: 3,
      color: '#84368d', // purple
    },
    {
      route: '/publications',
      title: 'Publications',
      navIndex: 4,
    },
    {
      route: '/faq',
      title: 'FAQ',
      navIndex: 5,
    },
    {
      route: '/about',
      title: 'About',
      // cardTitle: 'About',
      image: 'assets/images/gwas.svg',
      navIndex: 6,
    },
  ];

  // load data for publications tab
  const dispatch = useDispatch();
  const publicationsState = useSelector((state) => state.publications);
  const mergePublications = (state) => dispatch(actions.mergeState(state));

  useEffect(() => {
    const getData = async () => {
      const data = await (await fetch(`api/getPublications`)).json();

      const reducer = (acc, column) => [
        ...acc,
        {
          Header: column,
          accessor: column,
          id: column,
          Cell: (e) => {
            if (
              column == 'Title' &&
              e.row.values['DOI'] &&
              e.row.values['DOI'] != 'NA'
            ) {
              return (
                <a href={e.row.values['DOI']} target="_blank" rel="noreferrer">
                  {e.value}
                </a>
              );
            } else if (
              column == 'Name' &&
              e.row.values['Github'] &&
              e.row.values['Github'] != 'NA'
            ) {
              return (
                <a
                  href={e.row.values['Github']}
                  target="_blank"
                  rel="noreferrer"
                >
                  {e.value}
                </a>
              );
            } else {
              return e.value || '';
            }
          },
        },
      ];

      mergePublications({
        orA: {
          columns: [
            ...new Set(
              ...data['Original Research A'].map((row) => Object.keys(row))
            ),
          ].reduce(reducer, []),
          data: data['Original Research A'],
        },
        orB: {
          columns: [
            ...new Set(
              ...data['Orignal Research B'].map((row) => Object.keys(row))
            ),
          ].reduce(reducer, []),
          data: data['Orignal Research B'],
        },
        rp: {
          columns: [
            ...new Set(...data['Review Paper'].map((row) => Object.keys(row))),
          ].reduce(reducer, []),
          data: data['Review Paper'],
        },
        cm: {
          columns: [
            ...new Set(
              ...data['Computational Methods'].map((row) => Object.keys(row))
            ),
          ].reduce(reducer, []),
          data: data['Computational Methods'],
        },
      });
    };

    // get data on inital page load
    if (!publicationsState.orA.data) getData();
  }, [publicationsState]);

  return (
    <Router>
      <ErrorModal />
      <SuccessModal />
      <Navbar links={links} />
      <Route path="/" exact={true} render={(_) => <Home links={links} />} />
      <Route path="/about" component={About} />
      <Route path="/visualization/:type?/:id?" component={Visualization} />
      <Route path="/catalog" component={Catalog} />
      <Route path="/exploration/:exampleName?" component={Exposure} />
      <Route path="/refitting" component={Refitting} />
      <Route path="/association" component={Association} />
      <Route path="/publications" component={Publications} />
      <Route path="/faq" component={Faq} />
    </Router>
  );
}
