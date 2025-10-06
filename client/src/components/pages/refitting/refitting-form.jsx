import { useState, useEffect } from 'react';
import { Form, Button, Row, Col } from 'react-bootstrap';
import { useForm, Controller } from 'react-hook-form';
import { useHistory, useParams } from 'react-router-dom';
import SelectForm from '../../controls/select/selectHookForm';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { useSelector, useDispatch } from 'react-redux';

export default function RefittingForm() {
  const [submitted, setSubmitted] = useState(false);
  const [loadingUpload, setLoadingUpload] = useState(false);
  const [loadingSubmit, setLoadingSubmit] = useState(false);

  const defaultValues = {
    source: 'public',
    signatureType: { label: 'SBS', value: 'SBS' },
    referenceGenome: { label: 'hg38 (GRCh38)', value: 'hg38' },
    mafFile: null,
    genomicFile: null,
    clinicalFile: null,
    email: '',
    jobName: '',
  };

  const {
    control,
    register,
    handleSubmit,
    reset: resetForm,
    setValue,
    watch,
    formState: { errors },
  } = useForm({ defaultValues: defaultValues });

  const {
    source,
    signatureType,
    referenceGenome,
    mafFile,
    genomicFile,
    clinicalFile,
    email,
    jobName,
  } = watch();

  const signatureTypeOptions = [
    { label: 'SBS', value: 'SBS' },
    { label: 'DBS', value: 'DBS' },
  ];

  const referenceGenomeOptions = [
    { label: 'hg19 (GRCh37)', value: 'hg19' },
    { label: 'hg38 (GRCh38)', value: 'hg38' },
  ];

  const onSubmit = async (data) => {
    setLoadingSubmit(true);
    try {
      console.log('Refitting form submitted:', data);
      // Here you would implement the actual submission logic
      await new Promise(resolve => setTimeout(resolve, 2000)); // Simulate API call
      setSubmitted(true);
    } catch (error) {
      console.error('Submission error:', error);
    } finally {
      setLoadingSubmit(false);
    }
  };

  const handleReset = () => {
    resetForm(defaultValues);
    setSubmitted(false);
  };

  const handleFileUpload = async (file, fieldName) => {
    setLoadingUpload(true);
    try {
      // Simulate file upload
      await new Promise(resolve => setTimeout(resolve, 1000));
      setValue(fieldName, file);
    } catch (error) {
      console.error('File upload error:', error);
    } finally {
      setLoadingUpload(false);
    }
  };

  return (
    <div className="p-3 bg-white border rounded">
      <Form onSubmit={handleSubmit(onSubmit)}>
        <LoadingOverlay active={loadingUpload || loadingSubmit} />
        
        <div style={{ maxHeight: '900px', overflow: 'hidden auto' }}>
          <div className="border rounded p-2 mb-3">
            <Form.Group>
              <Form.Label className="mr-4">Data Source</Form.Label>
              <Controller
                name="source"
                control={control}
                render={({ field }) => (
                  <Form.Check
                    {...field}
                    inline
                    id="radioPublic"
                    type="radio"
                    label={<span className="font-weight-normal">Public</span>}
                    value={'public'}
                    checked={field.value === 'public'}
                    disabled={submitted}
                  />
                )}
              />
              <Controller
                name="source"
                control={control}
                render={({ field }) => (
                  <Form.Check
                    {...field}
                    inline
                    id="radioUser"
                    type="radio"
                    label={<span className="font-weight-normal">User</span>}
                    value={'user'}
                    checked={field.value === 'user'}
                    disabled={submitted}
                  />
                )}
              />
            </Form.Group>

            {source === 'public' ? (
              <div key="public">
                <p className="text-muted">
                  Use pre-loaded datasets available on the website for analysis.
                </p>
              </div>
            ) : (
              <div key="user">
                <SelectForm
                  className="mb-2"
                  name="signatureType"
                  label="Signature Type"
                  disabled={submitted}
                  options={signatureTypeOptions}
                  control={control}
                />
                
                <SelectForm
                  className="mb-2"
                  name="referenceGenome"
                  label="Reference Genome"
                  disabled={submitted}
                  options={referenceGenomeOptions}
                  control={control}
                />

                <Form.Group className="mb-2">
                  <Form.Label>Upload MAF File</Form.Label>
                  <Form.Control
                    type="file"
                    accept=".maf,.txt,.tsv"
                    disabled={submitted}
                    onChange={(e) => handleFileUpload(e.target.files[0], 'mafFile')}
                  />
                  <Form.Text className="text-muted">
                    Contains SBS or DBS information for samples
                  </Form.Text>
                </Form.Group>

                <Form.Group className="mb-2">
                  <Form.Label>UploadGenomic File</Form.Label>
                  <Form.Control
                    type="file"
                    accept=".bed,.txt,.tsv"
                    disabled={submitted}
                    onChange={(e) => handleFileUpload(e.target.files[0], 'genomicFile')}
                  />
                  <Form.Text className="text-muted">
                    Defines genomic regions targeted by sequencing panels
                  </Form.Text>
                </Form.Group>

                <Form.Group className="mb-2">
                  <Form.Label>UploadClinical File</Form.Label>
                  <Form.Control
                    type="file"
                    accept=".txt,.tsv,.csv"
                    disabled={submitted}
                    onChange={(e) => handleFileUpload(e.target.files[0], 'clinicalFile')}
                  />
                  <Form.Text className="text-muted">
                    Sample ID, sequencing panel ID, and cancer type
                  </Form.Text>
                </Form.Group>
              </div>
            )}
          </div>

          <div className="border rounded p-2 mb-3">
            

            <Form.Group className="mb-2">
              <Form.Label className="required">Email</Form.Label>
              <Form.Control
                {...register('email', { 
                  required: 'Email is required',
                  pattern: {
                    value: /^[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}$/i,
                    message: 'Invalid email address'
                  }
                })}
                type="email"
                placeholder="Enter your email address"
                disabled={submitted}
                isInvalid={!!errors.email}
              />
              <Form.Control.Feedback type="invalid">
                {errors.email?.message}
              </Form.Control.Feedback>
              <Form.Text className="text-muted">
                You will receive an email when the analysis is complete
              </Form.Text>
            </Form.Group>

            <Form.Group className="mb-2">
              <Form.Label>Job Name</Form.Label>
              <Form.Control
                {...register('jobName')}
                type="text"
                placeholder="Enter a descriptive name for your analysis"
                disabled={submitted}
              />
            </Form.Group>
          </div>

          <Row className="mt-3">
            <Col>
              <Button
                variant="secondary"
                onClick={handleReset}
                disabled={loadingSubmit}
                className="w-100"
              >
                Reset
              </Button>
            </Col>
            <Col>
              <Button
                variant="primary"
                type="submit"
                disabled={submitted || loadingSubmit}
                className="w-100"
              >
                {loadingSubmit ? 'Submitting...' : 'Submit'}
              </Button>
            </Col>
          </Row>
        </div>
      </Form>
    </div>
  );
}
