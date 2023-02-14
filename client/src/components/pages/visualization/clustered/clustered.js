import { useSelector } from 'react-redux';
import Description from '../../../controls/description/description';
import { NavHashLink } from 'react-router-hash-link';
import ClusteredPlot from './clustered-plot';
import ClusteredForm from './clustered-form';

export default function ClusteredIdentification() {
  const store = useSelector((state) => state.visualization);
  const { cluster, id, inputFormat } = store.userForm;

  return (
    <div className="bg-white border rounded" style={{ minHeight: '500px' }}>
      <Description
        className="p-3"
        less="This analysis identifies the kataegis events from a VCF file input."
        more={
          <>
            <span>
              Kataegis is a localized substitution hypermutation event, often
              characterized by clusters of C&#x3c;T and/or C&#x3c;G mutations,
              commonly at TpCpN trinucleotides (APOBEC mutations). Click{' '}
              <NavHashLink to="/faq#kataegis">here</NavHashLink> for additional
              information about kataegis.
            </span>
            <p className="mt-3">
              To identify kataegis, input the [Minimum Number of Mutations]
              required for kataegis, the [Maximum Distance] between one mutation
              and the next within a given cluster of mutations being considered
              for kataegis, and a [Chromosome] to be highlighted in the rainfall
              plot. By default, all chromosomes will be shown for the kataegis
              identification.
            </p>
          </>
        }
      />
      <hr />
      {cluster && id && inputFormat.value == 'vcf' && (
        <>
          <ClusteredForm />
          <hr />
          <ClusteredPlot />
        </>
      )}
      {inputFormat.value == 'vcf' && !cluster && (
        <div className="p-3">
          Calculation must be performed with the "Clustered Mutational Analysis"
          option enabled in the inital submission.
        </div>
      )}
      {inputFormat.value != 'vcf' && (
        <div className="p-3">Analysis is only available for VCF file input.</div>
      )}
    </div>
  );
}
