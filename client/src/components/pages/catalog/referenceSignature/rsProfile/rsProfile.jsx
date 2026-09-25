import RsProfileFormPlot from '@/components/pages/catalog/referenceSignature/rsProfile/rsProfile-form-plot';

export default function RsProfile() {
  return (
    <div>
      <div className="p-3">
        Enter any [Signature Source], [Profile Name], [Reference Signature Set],
        [Experimental Strategy], and [Signature Name] below to visualize the
        mutational signature profile.
      </div>
      <hr />
      <RsProfileFormPlot />
    </div>
  );
}
