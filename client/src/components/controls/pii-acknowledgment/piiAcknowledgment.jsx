import { Form } from 'react-bootstrap';
import { Controller } from 'react-hook-form';

// Checkbox the user must tick before submitting their own data. Gates the
// form's action button; wired through react-hook-form so Reset clears it.
export default function PiiAcknowledgment({
  id,
  control,
  name = 'piiAcknowledged',
  disabled = false,
  className = 'mb-3',
}) {
  return (
    <Controller
      name={name}
      control={control}
      render={({ field }) => (
        <Form.Group controlId={id} className={className}>
          <Form.Check
            type="checkbox"
            checked={!!field.value}
            onChange={(e) => field.onChange(e.target.checked)}
            onBlur={field.onBlur}
            disabled={disabled}
            label={
              <span className="font-weight-normal">
                I confirm that this data does not contain personally
                identifiable information (PII) or protected health information
                (PHI).
              </span>
            }
          />
        </Form.Group>
      )}
    />
  );
}
