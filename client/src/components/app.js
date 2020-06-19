import React from 'react';
import { HashRouter as Router, Route } from 'react-router-dom';
import Navbar from './controls/navbar/navbar';
import Home from './pages/home/home';
import About from './pages/about/about';
import Visualize from './pages/visualize/visualize';
import Explore from './pages/explore/explore';
import Refitting from './pages/refitting/refitting';
import Association from './pages/association/association';
// import { ErrorModal } from './controls/error-modal/error-modal';

export default function App() {
  const links = [
    {
      route: '/visualize',
      action: 'Visualize',
      title: 'Visualize',
      cardTitle: 'Signature Visualization',
      cardText: 'Visualize mutational profiles',
      image: 'assets/images/gwas.svg',
      navIndex: 0,
    },
    {
      route: '/explore',
      action: 'Explore',
      title: 'Explore',
      cardTitle: 'Signature Exploring',
      cardText: 'Signature Exploring',
      image: 'assets/images/phenotypes.svg',
      navIndex: 1,
    },
    {
      route: '/refitting',
      action: 'Refitting',
      title: 'Refitting',
      cardTitle: 'Signature Refitting',
      cardText: 'Signature Refitting',
      image: 'assets/images/downloads.svg',
      navIndex: 2,
    },
    {
      route: '/association',
      action: 'Association',
      title: 'Association',
      cardTitle: 'Signature Association',
      cardText: 'Signature Association',
      image: 'assets/images/downloads.svg',
      navIndex: 3,
    },
    {
      route: '/about',
      title: 'About',
      // cardTitle: 'About',
      image: 'assets/images/about.svg',
      navIndex: 4,
    },
  ];

  return (
    <Router>
      {/* <ErrorModal /> */}
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
