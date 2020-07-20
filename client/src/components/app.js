import React from 'react';
import { HashRouter as Router, Route } from 'react-router-dom';
import { Navbar } from './controls/navbar/navbar';
import Home from './pages/home/home';
import About from './pages/about/about';
import Visualize from './pages/visualize/visualize';
import Explore from './pages/explore/explore';
import Refitting from './pages/refitting/refitting';
import Association from './pages/association/association';
import { ErrorModal } from './controls/error-modal/error-modal';

export default function App() {
  const links = [
    {
      route: '/visualize',
      action: 'Visualize',
      title: 'Visualize',
      cardTitle: 'Signature Visualization',
      cardText: 'Visualize mutational profiles',
      description: 'Interactively and comprehensively visualize mutation signature in both sample and study level, including different type and level of mutational profiles (SBS/INDEL/DBS/SV/CNV), PCA components and different mutational feature (kataegis mutation, mutation quality, drive gene mutation etc).',
      image: 'assets/images/visualize.png',
      navIndex: 0,
    },
    {
      route: '/explore',
      action: 'Explore',
      title: 'Explore',
      cardTitle: 'Signature Exploring',
      cardText: 'Signature Exploring',
      description: 'Systematically explore any reference or update to date published signatures with different profiles, version and etiology (endogenous vs Exogenous). Intergratively explore the landscape of signature exposure in different genomic studies, including TCGA, PCAWG, and our Sherlock-Lung study.',
      image: 'assets/images/explore.png',
      navIndex: 1,
    },
    {
      route: '/refitting',
      action: 'Refitting',
      title: 'Refitting',
      cardTitle: 'Signature Refitting',
      cardText: 'Signature Refitting',
      description: 'Comprehensively evaluate the accuracy of mutational signature based on different statistical variables (Cosine similarity, BIC, L2 norms etc) and re-decompsite signatures using different algorithms (SigProfiler, deconstructsig, bootstrapping method).',
      image: 'assets/images/refitting.png',
      navIndex: 2,
    },
    {
      route: '/association',
      action: 'Association',
      title: 'Association',
      cardTitle: 'Signature Association',
      cardText: 'Signature Association',
      description: 'Systematically analyze and visualize the association between mutational signature exposure and genomic or epigenomic features or other sample based variables (such as clinical data) in different genomic studies.',
      image: 'assets/images/association.png',
      navIndex: 3,
    },
    {
      route: '/about',
      title: 'About',
      // cardTitle: 'About',
      image: 'assets/images/gwas.svg',
      navIndex: 4,
    },
  ];

  return (
    <Router>
      <ErrorModal />
      <Navbar links={links} />
      <Route path="/" exact={true} render={(_) => <Home links={links} />} />
      <Route path="/about" component={About} />
      <Route path="/visualize" component={Visualize} />
      <Route path="/explore" component={Explore} />
      <Route path="/refitting" component={Refitting} />
      <Route path="/association" component={Association} />
    </Router>
  );
}
