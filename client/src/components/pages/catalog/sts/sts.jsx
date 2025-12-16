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
  const mergeSTS = (state) => dispatch(actions.mergeCatalog({ sts: state }));

  const { etiology, signature } = useSelector((state) => state.catalog.sts);

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
      const sortedEtiologies = [...new Set(data.map((e) => e.etiology))].sort(
        etiologySort
      );
      mergeSTS({
        etiology: sortedEtiologies[0],
      });
    }
  }, [data, etiology, signature]);

  return (
    <Container fluid className="p-4 bg-white border rounded">
      <h4 className="text-center mb-3">
        Targeted Sequencing Signatures Catalog
      </h4>
      <p>
        This webpage presents a pan-cancer catalog of mutational signatures
        derived from targeted sequencing (TS) data using the SATS method. The
        repertoire plots below display the frequencies of single base
        substitution (SBS) and double base substitution (DBS) signatures
        identified across 111,711 tumors spanning 23 cancer types in AACR
        Project GENIE (version 13.0). You may move the computer cursor over each
        color block or dot to view its frequency.
      </p>
      <p>
        In addition to signature frequencies, we provide the proposed etiologies
        for each mutational signature along with its profile plot. These
        profiles may differ slightly from those obtained through whole-genome or
        whole-exome sequencing (WGS/WES). WGS/WES signatures are derived from
        tumor mutational counts (TMC), whereas the signatures shown here are
        based on tumor mutational burden (TMB) estimated from TS panels.
        Consequently, WGS/WES- and TS-based signature profiles are similar but
        not identical (e.g., SBS5).
      </p>

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
          {error?.data || 'An error occurred while retrieving STS signatures'}
        </div>
      )}
    </Container>
  );
}
