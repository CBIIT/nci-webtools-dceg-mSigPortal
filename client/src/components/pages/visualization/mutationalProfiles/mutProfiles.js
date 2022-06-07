import { Suspense } from "react";
import { Alert, Container } from "react-bootstrap";
import Loader from "../../../controls/loader/loader";
import ErrorBoundary from "../../../controls/errorBoundary/error-boundary";
import MutProfilesForm from "./mutProfiles-form";
import MutProfilesPlot from "./mutProfiles-plot";

export default function MutationalProfiles(props) {
  return (
    <Container
      fluid
      className="bg-white border rounded p-3 text-center"
      style={{ minHeight: 500 }}
      {...props}
    >
      <ErrorBoundary
        fallback={
          <Alert variant="danger">
            An internal error prevented plots from loading. Please contact the
            website administrator if this problem persists.
          </Alert>
        }
      >
        <Suspense fallback={<Loader message="Loading Plot" />}>
          <MutProfilesForm />
          <MutProfilesPlot />
        </Suspense>
      </ErrorBoundary>
    </Container>
  );
}