import { useEffect } from 'react';
import { Container } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../../services/store/catalog';
import STSEtiologyOptions from './stsEtiologyOptions';
import STSSignatureOptions from './stsSignatureOptions';
import STSSignatureInfo from './stsSignatureInfo';
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

  // automatically choose first etiology
  useEffect(() => {
    if (data && data[0]?.category == category && !etiology) {
      mergeSTS({
        etiology: [...new Set(data.map((e) => e.etiology))].sort()[0],
      });
    }
  }, [data, etiology, signature]);

  return (
    <Container fluid className="p-4 bg-white border rounded">
      <h5 className="separator mb-4">Signatures from Targeted Sequencing</h5>
      <p className="text-muted">
        Unlike COSMIC signatures, which are derived from tumor mutational counts (TMC) obtained from whole-genome or whole-exome sequencing (WGS/WES), signatures from targeted sequencing are derived from tumor mutational burden (TMB) estimated using targeted sequencing panels. Consequently, the profiles of TMC- and TMB-based signatures are similar but not identical (e.g., SBS5).
      </p>
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
