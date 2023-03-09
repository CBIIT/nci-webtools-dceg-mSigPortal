import Form from 'react-bootstrap/Form';
import Button from 'react-bootstrap/Button';
import InputGroup from 'react-bootstrap/InputGroup';

export default function HeaderSearch() {
  function handleSubmit(event) {
    event.preventDefault();

    const form = event.target;
    const site = window.location.host;
    const search = form.elements.q.value;
    const query = `site:${site} ${search}`;

    const searchParams = new URLSearchParams();
    searchParams.append('q', query);

    const searchUrl = `${form.action}?${searchParams}`;
    window.open(searchUrl, '_blank');
  }

  return (
    <Form
      className="d-flex align-items-stretch mb-4 mb-md-0"
      role="search"
      action="https://www.google.com/search"
      onSubmit={handleSubmit}
    >
      <InputGroup className="border-white">
        <Form.Control
          className="search-control"
          type="search"
          placeholder="search"
          aria-label="search"
          name="q"
        />
        <InputGroup.Append>
          <Button
            variant="outline-secondary"
            className="search-control-button"
            type="submit"
          >
            <i className="bi bi-search"></i>
            <span className="sr-only">submit</span>
          </Button>
        </InputGroup.Append>
      </InputGroup>
    </Form>
  );
}
