import React from 'react';
import { createRoot } from 'react-dom/client';
import App from './components/app';
import { Provider } from 'react-redux';
import { store } from './services/store';
import './index.scss';
import './ncids.scss';
import 'font-awesome/css/font-awesome.min.css';
import '@fontsource/montserrat';
import '@fontsource/inter';

createRoot(document.getElementById('root')).render(
  <Provider store={store}>
    <App />
  </Provider>
);
