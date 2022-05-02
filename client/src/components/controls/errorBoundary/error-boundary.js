import { Component } from 'react';

// Error boundaries currently have to be classes.
export default class ErrorBoundary extends Component {
  state = {
    hasError: false,
    error: null,
  };

  static getDerivedStateFromError(error) {
    return {
      hasError: true,
      error,
    };
  }

  render() {
    const { fallback, children } = this.props;
    const { hasError } = this.state;
    return hasError ? fallback : children;
  }
}
