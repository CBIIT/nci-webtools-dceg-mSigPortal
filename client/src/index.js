import React from 'react';
import ReactDOM from 'react-dom';
import App from './components/app';
import { Provider } from 'react-redux';
import { store } from './services/visualSlice';
import './index.scss';
import 'font-awesome/css/font-awesome.min.css';


ReactDOM.render(
    <Provider store={store}>
      <App />
    </Provider>, 
    document.getElementById('root')
);