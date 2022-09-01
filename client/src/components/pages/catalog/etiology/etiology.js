import { useEffect } from 'react';
import { Row, Col, Button } from 'react-bootstrap';
import { useSelector, useDispatch } from 'react-redux';
import axios from 'axios';
import { actions as catalogActions } from '../../../../services/store/catalog';
import { actions as modalActions } from '../../../../services/store/modal';
import { getJSON } from '../../../../services/utils';
import General from './general';
import CancerSpecific from './cancerSpecific';
import Therapies from './therapies';
import './etiology.scss';

const actions = { ...catalogActions, ...modalActions };

export default function Etiology() {
  const dispatch = useDispatch();
  const mergeEtiology = (state) =>
    dispatch(actions.mergeCatalog({ etiology: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const {
    category,
    selectedSignature,
    tissue,
    refSig,
    study,
    data,
    thumbnails,
    tissueThumbnails,
    refSigThumbnails,
    profileURL,
    exposureURL,
    strandbiasURL,
    tissueURL,
    refSigURL,
    loading,
  } = useSelector((state) => state.catalog.etiology);

  const categories = [
    {
      name: 'Cosmic Mutational Signatures (v3.2)',
      author: 'Alexandrov et al., 2021',
      etiologyTitle: 'Proposed Etiologies',
      file: 'Etiology_cosmic.json',
    },
    {
      name: 'Environmental Mutagenesis',
      author: 'Nik-Zainal et al., 2019',
      etiologyTitle: 'Proposed Mutagens',
      file: 'Etiology_enviromental_mutagenesis.json',
    },
    {
      name: 'DNA Repair Gene Edits',
      author: 'Nik-Zainal et al., 2018 and 2021',
      etiologyTitle: 'Genes',
      file: 'Etiology_gene_edits.json',
    },
    {
      name: 'Cancer Specific Signatures',
      author: 'Nik-Zainal et al., 2020',
      file: 'Etiology_cancer_specific_signatures.json',
    },
    {
      name: 'Cancer Therapies',
      author: 'Lopez-Bigas et al., 2019',
      etiologyTitle: 'Therapies',
      file: 'Etiology_cancer_therapies.json',
    },
    {
      name: 'Others',
      etiologyTitle: 'Proposed Etiologies',
      file: 'Etiology_others.json',
    },
  ];

  useEffect(() => {
    dispatch(actions.mergeCatalog({ catalog: { displayTab: 'etiology' } }));
  }, []);

  useEffect(() => {
    const getData = async () => {
      mergeEtiology({ loading: true });
      try {
        const file = categories.filter(({ name }) => name == category)[0].file;
        const data = await getJSON(`Etiology/${file}`);

        const etiologyOptions = [
          ...new Set(
            data.map(({ Etiology, Treatments }) => Etiology || Treatments)
          ),
        ];
        const studyOptions = [...new Set(data.map(({ Study }) => Study))];

        mergeEtiology({
          data: { [category]: data },
          etiology: etiologyOptions[0],
          study: studyOptions[0],
        });
      } catch (e) {
        mergeError(e.message);
        console.error(e);
      }
      mergeEtiology({ loading: false });
    };
    if (!data[category]) {
      if (!loading) getData();
    } else {
      const etiologyOptions = [
        ...new Set(
          data[category].map(
            ({ Etiology, Treatments }) => Etiology || Treatments
          )
        ),
      ];
      const studyOptions = [
        ...new Set(data[category].map(({ Study }) => Study)),
      ];

      mergeEtiology({
        etiology: etiologyOptions[0],
        study: studyOptions[0],
      });
    }
  }, [category, loading]);

  // check if plot exists
  useEffect(() => {
    const getImageS3 = (path) =>
      axios
        .post(
          'web/getImageS3',
          { path: `msigportal/Database/Etiology/${path}` },
          { responseType: 'blob' }
        )
        .then((response) => URL.createObjectURL(response.data))
        .catch((_) => '');

    const getPlots = async () => {
      if (profileURL) URL.revokeObjectURL(profileURL);
      if (exposureURL) URL.revokeObjectURL(exposureURL);
      if (strandbiasURL) URL.revokeObjectURL(strandbiasURL);

      const [sig, tmb, strandBias] = await [
        getImageS3(`Profile/${fixFile(selectedSignature)}.svg`),
        getImageS3(`Exposure/${fixFile(`${selectedSignature}_${study}`)}.svg`),
        getImageS3(`Profile_StrandBias/${fixFile(selectedSignature)}.svg`),
      ];
      mergeEtiology({
        profileURL: sig,
        exposureURL: tmb,
        strandbiasURL: strandBias,
      });
    };

    const getCancerSpecificPlots = async () => {
      if (tissueURL) URL.revokeObjectURL(tissueURL);
      if (refSigURL) URL.revokeObjectURL(refSigURL);
      if (exposureURL) URL.revokeObjectURL(exposureURL);
      const [tissuePlot, refSigPlot, exposurePlot] = await Promise.all([
        getImageS3(`Profile/${fixFile(tissue)}.svg`),
        getImageS3(`Profile/${fixFile(refSig)}.svg`),
        getImageS3(`Exposure/${fixFile(tissue)}_PCAWG.svg`),
      ]);
      mergeEtiology({
        tissueURL: tissuePlot,
        refSigURL: refSigPlot,
        exposureURL: exposurePlot,
      });
    };

    if (selectedSignature) getPlots();
    else if (tissue && refSig) getCancerSpecificPlots();
  }, [selectedSignature, study, tissue, refSig]);

  // get thumbnails for standard categories
  useEffect(() => {
    if (
      data[category] &&
      data[category].length &&
      !thumbnails[category] &&
      category != 'Cancer Specific Signatures' &&
      category != 'Cancer Therapies' &&
      !loading
    ) {
      const signatures = data[category]
        .slice()
        .filter(
          (value, index, array) =>
            array.findIndex(
              (t) => t['Signature Name'] === value['Signature Name']
            ) === index ||
            array.findIndex((t) => t['Signature'] === value['Signature']) ===
              index
        )
        .sort(naturalSort)
        .sort(profileSort);

      getThumbnails(signatures);
    }
  }, [data, category, loading]);

  // get tissue thumbnails
  useEffect(() => {
    if (
      data[category] &&
      data[category].length &&
      !tissueThumbnails.length &&
      category == 'Cancer Specific Signatures'
    ) {
      const uniqueTissues = Object.values(
        data[category].reduce((c, e) => {
          if (!c[e['Tissue Specific Signature']])
            c[e['Tissue Specific Signature']] = e;
          return c;
        }, {})
      ).sort(naturalSort);

      getTissueThumbnails(uniqueTissues);
    }
  }, [data, category]);

  // get refsig thumbnails
  useEffect(() => {
    if (
      data[category] &&
      data[category].length &&
      !refSigThumbnails.length &&
      category == 'Cancer Specific Signatures'
    ) {
      const refsig = data[category].slice().sort(naturalSort);

      getRefSigThumbnails(refsig);
    }
  }, [data]);

  // get therapy thumbnails
  useEffect(() => {
    if (
      data[category] &&
      data[category].length &&
      !thumbnails[category] &&
      category == 'Cancer Therapies' &&
      !loading
    ) {
      const signatures = data[category]
        .slice()
        .filter(
          (value, index, array) =>
            array.findIndex((t) => t['Signature'] === value['Signature']) ===
            index
        )
        .sort(naturalSort)
        .sort(profileSort);

      getTherapyThumbnails(signatures);
    }
  }, [data, category, loading]);

  function profileSort(a, b) {
    const sigOrder = [/SBS/g, /DBS/g, /ID/g];
    let c = 0,
      d = 0;

    sigOrder.forEach((profile, i) => {
      const sigA = a['Signature Name'] || a.Signature;
      const sigB = b['Signature Name'] || b.Signature;
      if (sigA.match(profile)) c = i;
      if (sigB.match(profile)) d = i;
    });

    return c - d;
  }

  // replace forward slash with colon for proper fs traversal
  function fixFile(filename) {
    return filename.replace(/\//g, ':');
  }

  function naturalSort(a, b) {
    const sigA = a['Signature Name'] || a.Signature || false;
    const sigB = b['Signature Name'] || b.Signature || false;
    if (sigA)
      return sigA.localeCompare(sigB, undefined, {
        numeric: true,
        sensitivity: 'base',
      });
    else if (a['Tissue Specific Signature'])
      return a['Tissue Specific Signature'].localeCompare(
        b['Tissue Specific Signature'],
        undefined,
        {
          numeric: true,
          sensitivity: 'base',
        }
      );
    else {
      return 0;
    }
  }

  function getCategories() {
    return (
      <Row className="justify-content-center mb-3">
        {categories.map(({ name, author, file }) => (
          <Col key={name} lg="2" md="3" sm="4" className="mb-3 d-flex">
            <Button
              size="sm"
              variant="dark"
              onClick={
                file && name != category
                  ? async () => {
                      mergeEtiology({
                        category: name,
                        etiology: '',
                        selectedSignature: '',
                        study: '',
                      });
                    }
                  : () => {}
              }
              className={category != name ? 'disabled' : ''}
              block
            >
              <div>
                <div>{name}</div>
                <div>{author}</div>
              </div>
            </Button>
          </Col>
        ))}
      </Row>
    );
  }

  async function getImageBatch(keyArr) {
    return axios
      .post('web/getImageS3Batch', { keys: keyArr })
      .then(({ data }) => data);
  }

  async function getThumbnails(signatures) {
    try {
      const keys = signatures.map(
        (signature) =>
          `msigportal/Database/Etiology/Profile_logo/${fixFile(
            signature['Signature Name'] || signature.Signature
          )}.svg`
      );

      const images = await getImageBatch(keys);

      const thumbnails = signatures.map((signature, i) => ({
        Etiology: signature.Etiology,
        Study: signature.Study,
        signatureName: signature['Signature Name'] || signature.Signature,
        thumbnailURL: images[i],
      }));

      mergeEtiology({ thumbnails: { [category]: thumbnails } });
    } catch (err) {
      mergeError(err.message);
      console.error(err);
    }
  }

  async function getTissueThumbnails(tissues) {
    try {
      const keys = tissues.map(
        (tissue) =>
          `msigportal/Database/Etiology/Profile_logo/${fixFile(
            tissue['Tissue Specific Signature']
          )}.svg`
      );

      const images = await getImageBatch(keys);

      const thumbnails = tissues.map((tissue, i) => ({
        Etiology: tissue.Etiology,
        'Tissue Specific Signature': tissue['Tissue Specific Signature'],
        'Ref Signature': tissue['Ref Signature'],
        thumbnailURL: images[i],
      }));

      mergeEtiology({ tissueThumbnails: thumbnails });
    } catch (err) {
      mergeError(err.message);
      console.error(err);
    }
  }

  async function getRefSigThumbnails(refSigs) {
    try {
      const keys = refSigs.map(
        (refSig) =>
          `msigportal/Database/Etiology/Profile_logo/${fixFile(
            refSig['Ref Signature']
          )}.svg`
      );

      const images = await getImageBatch(keys);

      const thumbnails = refSigs.map((refSig, i) => ({
        Etiology: refSig.Etiology,
        'Tissue Specific Signature': refSig['Tissue Specific Signature'],
        'Ref Signature': refSig['Ref Signature'],
        'RefSig Proportion': refSig['RefSig Proportion'],
        thumbnailURL: images[i],
      }));

      mergeEtiology({ refSigThumbnails: thumbnails });
    } catch (err) {
      mergeError(err.message);
      console.error(err);
    }
  }

  async function getTherapyThumbnails(signatures) {
    try {
      const keys = signatures.map(
        (signature) =>
          `msigportal/Database/Etiology/Profile_logo/${fixFile(
            signature.Signature
          )}.svg`
      );

      const images = await getImageBatch(keys);

      const thumbnails = signatures.map((signature, i) => ({
        Etiology: signature.Treatments,
        Study: signature.Study,
        signatureName: signature.Signature,
        thumbnailURL: images[i],
      }));

      mergeEtiology({ thumbnails: { [category]: thumbnails } });
    } catch (err) {
      mergeError(err.message);
      console.error(err);
    }
  }

  return (
    <div className="p-4 bg-white border rounded">
      <h5 className="separator">Categories</h5>
      {getCategories()}
      {category != 'Cancer Specific Signatures' &&
        category != 'Cancer Therapies' && <General categories={categories} />}

      {category == 'Cancer Specific Signatures' && <CancerSpecific />}
      {category == 'Cancer Therapies' && <Therapies />}
    </div>
  );
}
