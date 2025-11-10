import { useState, useEffect } from 'react';
import { Form, Button, Row, Col, Alert } from 'react-bootstrap';
import { useForm, Controller } from 'react-hook-form';
import { useHistory, useParams } from 'react-router-dom';
import SelectForm from '../../controls/select/selectHookForm';
import { useSelector, useDispatch } from 'react-redux';
import { useSubmitRefittingMutation } from './apiSlice';
import { actions as modalActions } from '../../../services/store/modal';
import { actions as refittingActions } from '../../../services/store/refitting';

export default function RefittingForm() {
  const { submitted, ...state } = useSelector((state) => state.refitting.main);
  const id = state.id || false;
  
  const [localSubmitted, setLocalSubmitted] = useState(false);
  const [jobId, setJobId] = useState(id);
  const [error, setError] = useState(null);
  const [success, setSuccess] = useState(null);

  // API hooks
  const [submitRefitting] = useSubmitRefittingMutation();

  const dispatch = useDispatch();
  const mergeState = (state) => dispatch(refittingActions.mergeRefitting({ main: state }));
  const mergeSuccess = (msg) =>
    dispatch(modalActions.mergeModal({ success: { visible: true, message: msg } }));

  // Computed state for form disabled status
  const isFormDisabled = submitted || localSubmitted;

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

  // Sync form values to Redux store so Instructions component can access them
  useEffect(() => {
    dispatch(refittingActions.mergeRefitting({ 
      userForm: { 
        signatureType,
        referenceGenome,
        email,
        jobName,
      }
    }));
  }, [signatureType, referenceGenome, email, jobName, dispatch]);

  const signatureTypeOptions = [
    { label: 'SBS', value: 'SBS' },
    { label: 'DBS', value: 'DBS' },
  ];

  const referenceGenomeOptions = [
    { label: 'hg19 (GRCh37)', value: 'hg19' },
    { label: 'hg38 (GRCh38)', value: 'hg38' },
  ];

  // Validate MAF file content matches signature type
  const validateMafFile = async (file, signatureType) => {
    if (!file) return { isValid: true };
    
    try {
      const text = await file.text();
      const lines = text.split('\n').filter(line => line.trim());
      
      if (lines.length < 2) {
        return { isValid: false, error: 'MAF file appears to be empty or invalid' };
      }
      
      // Check header contains required columns
      const header = lines[0].toLowerCase();
      const requiredColumns = ['variant_type', 'chromosome', 'start_position', 'reference_allele', 'tumor_seq_allele2', 'tumor_sample_barcode'];
      const missingColumns = requiredColumns.filter(col => !header.includes(col));
      
      if (missingColumns.length > 0) {
        return { 
          isValid: false, 
          error: `MAF file missing required columns: ${missingColumns.join(', ')}` 
        };
      }
      
      // Check variant types in data rows
      const dataLines = lines.slice(1, Math.min(lines.length, 100)); // Check first 100 data lines
      const variantTypeIndex = lines[0].split('\t').findIndex(col => 
        col.toLowerCase().includes('variant_type')
      );
      
      if (variantTypeIndex === -1) return { isValid: true }; // Skip validation if column not found
      
      const variantTypes = dataLines
        .map(line => line.split('\t')[variantTypeIndex])
        .filter(vt => vt && vt.trim())
        .map(vt => vt.trim().toUpperCase());
      
      const expectedType = signatureType === 'SBS' ? 'SNP' : 'DNP';
      const hasExpectedType = variantTypes.some(vt => vt === expectedType);
      const hasWrongType = variantTypes.some(vt => vt === (signatureType === 'SBS' ? 'DNP' : 'SNP'));
      
      if (hasWrongType && !hasExpectedType) {
        return {
          isValid: false,
          error: `This MAF file contains ${signatureType === 'SBS' ? 'DBS (DNP)' : 'SBS (SNP)'} variants but you selected ${signatureType} analysis. Please select the correct signature type or upload the appropriate MAF file.`
        };
      }
      
      return { isValid: true };
    } catch (error) {
      return { 
        isValid: false, 
        error: 'Error reading MAF file. Please check the file format.' 
      };
    }
  };

  const onSubmit = async (data) => {
    // Validate required files and fields
    if (!data.mafFile || !data.genomicFile || !data.clinicalFile) {
      setError('Please upload all required files (MAF, Genomic, and Clinical files)');
      return;
    }
    
    if (!data.jobName || data.jobName.trim() === '') {
      setError('Job name is required');
      return;
    }

    // Validate MAF file content matches signature type
    const mafValidation = await validateMafFile(data.mafFile, data.signatureType);
    if (!mafValidation.isValid) {
      setError(mafValidation.error);
      return;
    }

    // Clear any existing errors
    setError(null);
    setSuccess(null);

    // Generate a UUID for tracking (matching extraction format)
    const jobId = crypto.randomUUID();
    
    // Update state and navigate to status tab immediately (like extraction)
    setJobId(jobId);
    setLocalSubmitted(true);
    mergeState({ id: jobId, displayTab: 'status' });
    
    // Show success modal immediately (like extraction)
    mergeSuccess(
      'Most Jobs take a long time, you will receive an email when the refitting job is complete. It is safe to close the window now'
    );

    // Store job in localStorage (like extraction)
    const jobs = JSON.parse(localStorage.getItem('refitting-jobs') || '[]');
    // const newJob = {
    //   id: jobId,
    //   jobName: data.jobName,
    //   status: 'SUBMITTED',
    //   submittedAt: new Date().toISOString(),
    //   email: data.email
    // };

    const updatedJobs = [...jobs, jobId];
    localStorage.setItem('refitting-jobs', JSON.stringify(updatedJobs));
    console.log('Stored refitting job:', jobId);
    console.log('All refitting jobs:', updatedJobs);

    // Create FormData for file upload
    const formData = new FormData();
    formData.append('mafFile', data.mafFile);
    formData.append('genomicFile', data.genomicFile);
    formData.append('clinicalFile', data.clinicalFile);
    formData.append('genome', data.referenceGenome);
    formData.append('signatureType', data.signatureType);
    formData.append('email', data.email);
    formData.append('jobName', data.jobName);

    // Submit in background - errors will be handled by status component
    try {
      console.log('About to submit refitting job...');
      console.log('Job ID:', jobId);
      console.log('FormData contents:');
      for (let pair of formData.entries()) {
        console.log(`  ${pair[0]}:`, pair[1]);
      }
      
      const result = await submitRefitting({ id: jobId, formData }).unwrap();
      console.log('Refitting job submitted successfully:', result);
    } catch (error) {
      console.error('Submission error details:', error);
      console.error('Error status:', error.status);
      console.error('Error data:', error.data);
      console.error('Error message:', error.message);
      console.error('Full error object:', JSON.stringify(error, null, 2));
      // Don't show errors here - let status component handle them
    }
  };

  const handleReset = () => {
    resetForm(defaultValues);
    setLocalSubmitted(false);
    setJobId(null);
    setError(null);
    setSuccess(null);
    mergeState({ id: null, displayTab: 'instructions', submitted: false });
  };

  return (
    <div className="p-3 bg-white border rounded">
      <Form onSubmit={handleSubmit(onSubmit)}>
        
        {/* Show validation errors only, not submission errors */}
        {error && !localSubmitted && (
          <Alert variant="danger" dismissible onClose={() => setError(null)}>
            {error}
          </Alert>
        )}

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
                        disabled={isFormDisabled}
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
                        disabled={isFormDisabled}
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
              <Form.Label>Upload {signatureType} MAF File <span style={{ color: 'crimson' }}>*</span></Form.Label>
              <Controller
                name="mafFile"
                control={control}
                rules={{ required: `${signatureType} MAF file is required` }}
                render={({ field }) => (
                  <Form.File
                    {...field}
                    value={''} // set dummy value for file input
                    disabled={isFormDisabled}
                    id="mafFile"
                    label={
                      mafFile?.name || `Upload ${signatureType} MAF File...`
                    }
                    isInvalid={errors.mafFile}
                    feedback={`Please upload a ${signatureType} MAF file`}
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
                disabled={isFormDisabled}
                onClick={async () => {
                  const file = `${signatureType}_MAF_two_samples.txt`;
                  const path = import.meta.env.BASE_URL + '/assets/examples/refitting/' + file;
                  setValue(
                    'mafFile',
                    new File([await (await fetch(path)).blob()], file)
                  );
                }}
                className='p-0'
              >
                Load {signatureType} Example
              </Button>
              <Button
                variant="link"
                disabled={isFormDisabled}
                onClick={() => {
                  const file = `${signatureType}_MAF_two_samples.txt`;
                  const path = import.meta.env.BASE_URL + '/assets/examples/refitting/' + file;
                  const link = document.createElement('a');
                  link.href = path;
                  link.download = file;
                  link.click();
                }}
                className='p-0'
              >
                Download {signatureType} MAF Example
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
                    disabled={isFormDisabled}
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
                disabled={isFormDisabled}
                onClick={async () => {
                  const file = 'Genomic_information_sample.txt';
                  const path = import.meta.env.BASE_URL + '/assets/examples/refitting/' + file;
                  setValue(
                    'genomicFile',
                    new File([await (await fetch(path)).blob()], file)
                  );
                }}
                className='p-0'
              >
                Load Example
              </Button>
              <Button
                variant="link"
                disabled={isFormDisabled}
                onClick={() => {
                  const file = 'Genomic_information_sample.txt';
                  const path = import.meta.env.BASE_URL + '/assets/examples/refitting/' + file;
                  const link = document.createElement('a');
                  link.href = path;
                  link.download = file;
                  link.click();
                }}
                className='p-0'
              >
                Download Genomic Example
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
                    disabled={isFormDisabled}
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
                disabled={isFormDisabled}
                onClick={async () => {
                  const file = 'Clinical_sample.txt';
                  const path = import.meta.env.BASE_URL + '/assets/examples/refitting/' + file;
                  setValue(
                    'clinicalFile',
                    new File([await (await fetch(path)).blob()], file)
                  );
                }}
                className='p-0'
              >
                Load Example
              </Button>
              <Button
                variant="link"
                disabled={isFormDisabled}
                onClick={() => {
                  const file = 'Clinical_sample.txt';
                  const path = import.meta.env.BASE_URL + '/assets/examples/refitting/' + file;
                  const link = document.createElement('a');
                  link.href = path;
                  link.download = file;
                  link.click();
                }}
                className='p-0'
              >
                Download Clinical Example
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
                disabled={isFormDisabled}
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
              <Form.Label className="required">Job Name</Form.Label>
              <Form.Control
                {...register('jobName', { 
                  required: 'Job name is required' 
                })}
                type="text"
                placeholder="Enter a descriptive name for your analysis"
                disabled={isFormDisabled}
                isInvalid={!!errors.jobName}
              />
              <Form.Control.Feedback type="invalid">
                {errors.jobName?.message}
              </Form.Control.Feedback>
            </Form.Group>
          </div>

          <Row className="mt-3">
            <Col>
              <Button
                variant="secondary"
                onClick={handleReset}
                className="w-100"
              >
                Reset
              </Button>
            </Col>
            <Col>
              <Button
                variant="primary"
                type="submit"
                disabled={isFormDisabled}
                className="w-100"
              >
                Submit
              </Button>
            </Col>
          </Row>
        </div>
      </Form>
    </div>
  );
}
