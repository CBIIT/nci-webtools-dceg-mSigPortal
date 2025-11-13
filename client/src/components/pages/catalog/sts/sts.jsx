import { useEffect } from 'react';
import { Container } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../../services/store/catalog';
import STSEtiologyOptions from './stsEtiologyOptions';
import STSSignatureOptions from './stsSignatureOptions';
import STSSignatureInfo from './stsSignatureInfo';
import STSTissueDistribution from './STSTissueDistribution';
import { useEtiologyOptionsQuery } from '../etiology/apiSlice';
import '../etiology/etiology.scss';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

const actions = { ...catalogActions };

export default function STS() {
  const dispatch = useDispatch();
  const mergeSTS = (state) =>
    dispatch(actions.mergeCatalog({ sts: state }));

  const { etiology, signature } = useSelector(
    (state) => state.catalog.sts
  );

  // Fixed category to STS
  const category = 'STS';

  const { data, error, isFetching } = useEtiologyOptionsQuery({
    category: category,
  });

  // Custom sort function matching the one in STSEtiologyOptions
  function etiologySort(a, b) {
    const cnRegex = /^CN:/;
    const unkRegex = /^Unknown/;
    const unkOrder = [
      'Unknown chemotherapy treatment',
      'Unknown (clock-like signature)',
      'Unknown',
    ];

    if (a.match(cnRegex)) {
      return 1;
    } else if (b.match(cnRegex)) {
      return -1;
    } else if (a.match(unkRegex) || b.match(unkRegex)) {
      return unkOrder.indexOf(a) - unkOrder.indexOf(b);
    } else if (a)
      return a.localeCompare(b, undefined, {
        numeric: true,
        sensitivity: 'base',
      });
    else {
      return 0;
    }
  }

  // automatically choose first etiology (using the same sort as display)
  useEffect(() => {
    if (data && data[0]?.category == category && !etiology) {
      const sortedEtiologies = [...new Set(data.map((e) => e.etiology))].sort(etiologySort);
      mergeSTS({
        etiology: sortedEtiologies[0],
      });
    }
  }, [data, etiology, signature]);

  return (
    <Container fluid className="p-4 bg-white border rounded">
      <div>
        Unlike COSMIC signatures, which are derived from tumor mutational counts (TMC) obtained from whole-genome or whole-exome sequencing (WGS/WES), signatures from targeted sequencing are derived from tumor mutational burden (TMB) estimated using targeted sequencing panels. Consequently, the profiles of TMC- and TMB-based signatures are similar but not identical (e.g., SBS5).
      </div>

      {/* Tissue Distribution with both SBS and DBS plots */}
      <STSTissueDistribution selectedSignature={signature} />
      
      {data && data[0]?.category == category ? (
        <>
          <STSEtiologyOptions data={data} />
          <STSSignatureOptions data={data} />
          <STSSignatureInfo data={data} />
        </>
      ) : (
        <div>
          <LoadingOverlay active={isFetching} />
        </div>
      )}
      {error && (
        <div>
          {error?.data || 'An error occured while retrieving STS signatures'}
        </div>
      )}
    </Container>
  );
}
