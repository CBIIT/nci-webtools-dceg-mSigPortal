import React from 'react';
import { HashRouter as Router, Route } from 'react-router-dom';
import { Navbar } from './controls/navbar/navbar';
import Home from './pages/home/home';
import About from './pages/about/about';
import Visualization from './pages/visualization/visualization';
import Exploring from './pages/exploring/exploring';
import Refitting from './pages/refitting/refitting';
import Association from './pages/association/association';
import Publications from './pages/publications/publications';
import Faq from './pages/faq/faq';
import { ErrorModal } from './controls/error-modal/error-modal';
import { SuccessModal } from './controls/success-modal/success-modal';

export default function App() {
  const links = [
    {
      route: '/visualization',
      action: 'Visualization',
      title: 'Visualization',
      cardTitle: 'Signature Visualization',
      cardText: 'Visualize mutational profiles',
      description:
        'Interactively and comprehensively visualize mutation signature in both sample and study level, including different type and level of mutational profiles (SBS/INDEL/DBS/SV/CNV), PCA components and different mutational feature (kataegis mutation, mutation quality, drive gene mutation etc).',
      image: 'assets/images/visualize.png',
      navIndex: 0,
      color: '#fc8701', // orange
      examples: [],
      // examples: [{ title: 'VCF Example', folder: 'vcfExample' }],
    },
    {
      route: '/exploring',
      action: 'Exploring',
      title: 'Exploring',
      cardTitle: 'Signature Exploring',
      cardText: 'Signature Exploring',
      description:
        'Systematically explore any reference or update to date published signatures with different profiles, version and etiology (endogenous vs Exogenous). Intergratively explore the landscape of signature exposure in different genomic studies, including TCGA, PCAWG, and our Sherlock-Lung study.',
      image: 'assets/images/explore.png',
      navIndex: 1,
      color: '#2c71dd', // blue
      examples: [],
      dropdown: [
        { name: 'Etiology', path: 'etiology' },
        { name: 'Signature', path: 'signature' },
        { name: 'Exposure', path: 'exposure' },
        { name: 'Download', path: 'download' },
      ],
    },
    {
      route: '/refitting',
      action: 'Refitting',
      title: 'Refitting',
      cardTitle: 'Signature Refitting',
      cardText: 'Signature Refitting',
      description:
        'Comprehensively evaluate the accuracy of mutational signature based on different statistical variables (Cosine similarity, BIC, L2 norms etc) and re-decompsite signatures using different algorithms (SigProfiler, deconstructsig, bootstrapping method).',
      image: 'assets/images/refitting.png',
      navIndex: 2,
      color: '#689f39', // green
      examples: [],
    },
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
      examples: [],
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

  return (
    <Router>
      <ErrorModal />
      <SuccessModal />
      <Navbar links={links} />
      <Route path="/" exact={true} render={(_) => <Home links={links} />} />
      <Route path="/about" component={About} />
      <Route path="/visualization/:type?/:id?" component={Visualization} />
      <Route path="/exploring" component={Exploring} />
      <Route path="/refitting" component={Refitting} />
      <Route path="/association" component={Association} />
      <Route path="/publications" component={Publications} />
      <Route path="/faq" component={Faq} />
    </Router>
  );
}
