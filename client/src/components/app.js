import React, { useEffect } from 'react';
import { RecoilRoot } from 'recoil';
import { HashRouter as Router, Route } from 'react-router-dom';
import { NavbarCustom } from './controls/navbar/navbar';
import Home from './pages/home/home';
import About from './pages/about/about';
import Visualization from './pages/visualization/visualization';
import Catalog from './pages/catalog/catalog';
import Exposure from './pages/exposure/exposure';
import Refitting from './pages/refitting/refitting';
import Association from './pages/association/association';
import Extraction from './pages/extraction/extraction';
import Publications from './pages/publications/publications';
import Faq from './pages/faq/faq';
import APIAccess from './pages/apiAccess/apiAccess';
import { Header } from './controls/header/header';
import { Footer } from './controls/footer/footer';
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
      cardId: 'Catalog',
      name: 'catalog',
      cardTitle: 'Signature Catalog',
      cardText: 'Signature Catalog',
      description:
        'All existing human and mouse signatures based on different genome builds and algorithm versions.',
      image: 'assets/images/Catalog-Icon.svg',
      navIndex: 0,
      color: '#4b833a',
      showHomepage: true,
      about: '/about#catalog',
    },

    {
      route: '/visualization',
      action: 'Visualization',
      title: 'Visualization',
      cardId: 'Visualization',
      name: 'visualization',
      cardTitle: 'Signature Visualization',
      cardText: 'Visualize mutational profiles',
      description:
        'Allow identication of signature features at sample level and sicovery of new signatures.',
      image: 'assets/images/Visualization-Icon.svg',
      navIndex: 1,
      color: '#2c5b4e',
      showHomepage: true,
      about: '/about#visualization',
    },
    {
      route: '/extraction',
      action: 'Extraction',
      title: 'Extraction',
      cardId: 'Extraction',
      name: 'extraction',
      cardTitle: 'Signature Extraction',
      cardText: 'Extraction mutational profiles',
      description:
        'Extract and compare muational signatures using state-of-the-art algorithms.',
      image: 'assets/images/Extraction-Icon.svg',
      navIndex: 2,
      color: '#2f4a64',
      showHomepage: true,
      about: '/about#extraction',
    },
    {
      route: '/exploration',
      action: 'Exploration',
      title: 'Exploration',
      cardId: 'Exploration',
      name: 'exploration',
      cardTitle: 'Signature Exploration',
      cardText: 'Signature Exploration',
      description:
        'Explore etiological factors associated with signature at sample levels.    ',
      image: 'assets/images/Exploration-Icon.svg',
      navIndex: 3,
      color: '#5a4e2e',
      showHomepage: true,
      about: '/about#exploration',
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
      cardId: 'Association',
      name: 'association',
      cardTitle: 'Signature Association',
      cardText: 'Signature Association',
      description:
        'Analyze signature association with other genomic features and clinical data.',
      image: 'assets/images/Association-Icon.svg',
      navIndex: 4,
      color: '#7f282f',
      showHomepage: true,
      about: '/about#association',
    },
    {
      route: '/apiaccess',
      action: 'API Access',
      title: 'API Access',
      cardId: 'APIAccess',
      name: 'api',
      cardTitle: 'API Access',
      cardText: 'API Access',
      description:
        'Provide REST API for programmatic access to the mSigPortal data.',
      image: 'assets/images/API-Icon.svg',
      navIndex: 5,
      color: '#84368d',
      showHomepage: true,
      about: '/about#api',
    },
    {
      route: '/publications',
      title: 'Publications',
      navIndex: 6,
      showHomepage: false,
    },
    {
      route: '/faq',
      title: 'FAQ',
      navIndex: 7,
      showHomepage: false,
    },
    {
      route: '/about',
      title: 'About',
      // cardTitle: 'About',
      image: 'assets/images/gwas.svg',
      navIndex: 8,
      showHomepage: false,
    },
  ];

  // load data for publications tab
  const dispatch = useDispatch();
  const publicationsState = useSelector((state) => state.publications);
  const mergePublications = (state) => dispatch(actions.mergeState(state));

  useEffect(() => {
    const getData = async () => {
      const data = await (await fetch(`web/getPublications`)).json();

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
    <>
      <Header />
      <RecoilRoot>
        <Router>
          <ErrorModal />
          <SuccessModal />
          <NavbarCustom links={links} />
          <Route path="/" exact={true} render={(_) => <Home links={links} />} />
          <Route path="/about" component={About} />
          <Route path="/visualization/:type?/:id?" component={Visualization} />
          <Route path="/catalog" component={Catalog} />
          <Route path="/exploration/:exampleName?" component={Exposure} />
          <Route path="/refitting" component={Refitting} />
          <Route path="/association" component={Association} />
          <Route path="/extraction" component={Extraction} />
          <Route path="/publications" component={Publications} />
          <Route path="/faq" component={Faq} />
          <Route path="/apiaccess" component={APIAccess} />
        </Router>
      </RecoilRoot>
      <Footer />
    </>
  );
}
