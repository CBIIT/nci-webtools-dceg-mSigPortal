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
    signatureType: 'SBS',
    referenceGenome: 'hg19',
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

  return (
    <div className="p-3 bg-white border rounded">
      <Form onSubmit={handleSubmit(onSubmit)}>
        <LoadingOverlay active={loadingUpload || loadingSubmit} />
        
        <div style={{ maxHeight: '900px', overflow: 'hidden auto' }}>
          <div className="border rounded p-2 mb-3">
            <Form.Group className="mb-3">
              <Form.Label className="required mr-4">Signature Type</Form.Label>
              <Controller
                name="signatureType"
                control={control}
                rules={{ required: 'Signature type is required' }}
                render={({ field }) => (
                  <div>
                    {signatureTypeOptions.map((option) => (
                      <Form.Check
                        key={option.value}
                        type="radio"
                        id={`signatureType-${option.value}`}
                        label={<span className="font-weight-normal">{option.label}</span>}
                        value={option.value}
                        checked={field.value === option.value}
                        onChange={(e) => field.onChange(e.target.value)}
                        disabled={submitted}
                        inline
                      />
                    ))}
                  </div>
                )}
              />
              {errors.signatureType && (
                <Form.Text className="text-danger">
                  {errors.signatureType.message}
                </Form.Text>
              )}
            </Form.Group>
            
            <Form.Group className="mb-3">
              <Form.Label className="required mr-4">Reference Genome</Form.Label>
              <Controller
                name="referenceGenome"
                control={control}
                rules={{ required: 'Reference genome is required' }}
                render={({ field }) => (
                  <div>
                    {referenceGenomeOptions.map((option) => (
                      <Form.Check
                        key={option.value}
                        type="radio"
                        id={`referenceGenome-${option.value}`}
                        label={<span className="font-weight-normal">{option.label}</span>}
                        value={option.value}
                        checked={field.value === option.value}
                        onChange={(e) => field.onChange(e.target.value)}
                        disabled={submitted}
                        inline
                      />
                    ))}
                  </div>
                )}
              />
              {errors.referenceGenome && (
                <Form.Text className="text-danger">
                  {errors.referenceGenome.message}
                </Form.Text>
              )}
            </Form.Group>

            <Form.Group className="mb-2">
              <Form.Label>Upload MAF File <span style={{ color: 'crimson' }}>*</span></Form.Label>
              <Controller
                name="mafFile"
                control={control}
                rules={{ required: 'MAF file is required' }}
                render={({ field }) => (
                  <Form.File
                    {...field}
                    value={''} // set dummy value for file input
                    disabled={submitted}
                    id="mafFile"
                    label={
                      mafFile?.name || 'Upload MAF File...'
                    }
                    isInvalid={errors.mafFile}
                    feedback="Please upload a MAF file"
                    onChange={(e) => {
                      if (e.target.files.length) {
                        setValue('mafFile', e.target.files[0]);
                      }
                    }}
                    custom
                  />
                )}
              />
              <Button
                variant="link"
                disabled={submitted}
                onClick={async () => {
                  const file = 'SBS_MAF_two_samples.txt';
                  const path = 'assets/examples/refitting/' + file;
                  setValue(
                    'mafFile',
                    new File([await (await fetch(path)).blob()], file)
                  );
                }}
              >
                Load Example MAF File
              </Button>
              <Button
                variant="link"
                disabled={submitted}
                onClick={() => {
                  const file = 'SBS_MAF_two_samples.txt';
                  const path = 'assets/examples/refitting/' + file;
                  const link = document.createElement('a');
                  link.href = path;
                  link.download = file;
                  link.click();
                }}
              >
                Download Example MAF File
              </Button>
             
            </Form.Group>

            <Form.Group className="mb-2">
              <Form.Label>Upload Genomic File <span style={{ color: 'crimson' }}>*</span></Form.Label>
              <Controller
                name="genomicFile"
                control={control}
                rules={{ required: 'Genomic file is required' }}
                render={({ field }) => (
                  <Form.File
                    {...field}
                    value={''} // set dummy value for file input
                    disabled={submitted}
                    id="genomicFile"
                    label={
                      genomicFile?.name || 'Upload Genomic File...'
                    }
                    isInvalid={errors.genomicFile}
                    feedback="Please upload a genomic file"
                    onChange={(e) => {
                      if (e.target.files.length) {
                        setValue('genomicFile', e.target.files[0]);
                      }
                    }}
                    custom
                  />
                )}
              />
              <Button
                variant="link"
                disabled={submitted}
                onClick={async () => {
                  const file = 'Genomic_information_sample.txt';
                  const path = 'assets/examples/refitting/' + file;
                  setValue(
                    'genomicFile',
                    new File([await (await fetch(path)).blob()], file)
                  );
                }}
              >
                Load Example Genomic File
              </Button>
              <Button
                variant="link"
                disabled={submitted}
                onClick={() => {
                  const file = 'Genomic_information_sample.txt';
                  const path = 'assets/examples/refitting/' + file;
                  const link = document.createElement('a');
                  link.href = path;
                  link.download = file;
                  link.click();
                }}
              >
                Download Example Genomic File
              </Button>
              
            </Form.Group>

            <Form.Group className="mb-2">
              <Form.Label>Upload Clinical File <span style={{ color: 'crimson' }}>*</span></Form.Label>
              <Controller
                name="clinicalFile"
                control={control}
                rules={{ required: 'Clinical file is required' }}
                render={({ field }) => (
                  <Form.File
                    {...field}
                    value={''} // set dummy value for file input
                    disabled={submitted}
                    id="clinicalFile"
                    label={
                      clinicalFile?.name || 'Upload Clinical File...'
                    }
                    isInvalid={errors.clinicalFile}
                    feedback="Please upload a clinical file"
                    onChange={(e) => {
                      if (e.target.files.length) {
                        setValue('clinicalFile', e.target.files[0]);
                      }
                    }}
                    custom
                  />
                )}
              />
              <Button
                variant="link"
                disabled={submitted}
                onClick={async () => {
                  const file = 'Clinical_sample.txt';
                  const path = 'assets/examples/refitting/' + file;
                  setValue(
                    'clinicalFile',
                    new File([await (await fetch(path)).blob()], file)
                  );
                }}
              >
                Load Example Clinical File
              </Button>
              <Button
                variant="link"
                disabled={submitted}
                onClick={() => {
                  const file = 'Clinical_sample.txt';
                  const path = 'assets/examples/refitting/' + file;
                  const link = document.createElement('a');
                  link.href = path;
                  link.download = file;
                  link.click();
                }}
              >
                Download Example Clinical File
              </Button>
              
            </Form.Group>
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
