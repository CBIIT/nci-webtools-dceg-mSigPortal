import Spinner from 'react-bootstrap/Spinner';

export default function Loader({ message }) {
  return (
    <div className="loader">
      <Spinner variant="primary" animation="border" role="status" />
      <div>{message || 'Loading'}</div>
    </div>
  );
}
