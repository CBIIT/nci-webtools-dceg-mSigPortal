import { useEffect } from 'react';
import { Container } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import { actions as catalogActions } from '../../../../services/store/catalog';
import CategoryOptions from './categoryOptions';
import EtiologyOptions from './etiologyOptions';
import SignatureOptions from './signatureOptions';
import SignatureInfo from './signatureInfo';
import { useEtiologyOptionsQuery } from './apiSlice';
import './etiology.scss';
import { LoadingOverlay } from '../../../controls/loading-overlay/loading-overlay';

const actions = { ...catalogActions };

export default function Etiology() {
  const dispatch = useDispatch();
  const mergeEtiology = (state) =>
    dispatch(actions.mergeCatalog({ etiology: state }));

  const { category, etiology, signature } = useSelector(
    (state) => state.catalog.etiology
  );

  const categories = [
    {
      category: 'Cosmic',
      name: 'Cosmic Mutational Signatures (v3.3)',
      author: 'Alexandrov et al., 2021',
      etiologyTitle: 'Proposed Etiologies',
    },
    {
      category: 'CancerSpecificSignatures_2022',
      name: 'Cancer Reference Signatures',
      author: 'Nik-Zainal et al., 2022',
    },
    {
      category: 'EnviromentalMutagenesis',
      name: 'Environmental Mutagenesis',
      author: 'Nik-Zainal et al., 2019',
      etiologyTitle: 'Proposed Mutagens',
    },
    {
      category: 'GeneEdits',
      name: 'DNA Repair Gene Edits',
      author: 'Nik-Zainal et al., 2018 and 2021',
      etiologyTitle: 'Genes',
    },
    {
      category: 'CancerSpecificSignatures',
      name: 'Cancer Specific Signatures',
      author: 'Nik-Zainal et al., 2020',
    },
    {
      category: 'CancerTherapies',
      name: 'Cancer Therapies',
      author: 'Lopez-Bigas et al., 2019',
      etiologyTitle: 'Therapies',
    },
    {
      category: 'Others',
      name: 'Others',
      etiologyTitle: 'Proposed Etiologies',
    },
  ];

  const { data, error, isFetching } = useEtiologyOptionsQuery({
    category: category || categories[0].category,
  });

  // automatically choose first etiology
  useEffect(() => {
    if (data && data[0].category == category && !etiology) {
      mergeEtiology({
        etiology: [...new Set(data.map((e) => e.etiology))].sort()[0],
      });
    }
  }, [data, category, etiology, signature]);

  return (
    <Container fluid className="p-4 bg-white border rounded">
      <h5 className="separator">Categories</h5>
      <CategoryOptions categories={categories} />
      {data && data[0]?.category == category ? (
        <>
          <EtiologyOptions data={data} />
          <SignatureOptions data={data} />
          <SignatureInfo data={data} />
        </>
      ) : (
        <div>
          <LoadingOverlay active={isFetching} />
        </div>
      )}
      {error && (
        <div>
          {error?.data || 'An error occured while retrieving etiology options'}
        </div>
      )}
    </Container>
  );
}
