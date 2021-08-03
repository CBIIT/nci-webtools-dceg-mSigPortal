import React, { useEffect } from 'react';
import { Link } from 'react-router-dom';
import ProfilerSummary from './profilerSummary';
import MutationalProfiles from './mutationalProfiles';
import CosineSimilarity from './cosineSimilarity';
import MutationalPattern from './mutationalPattern';
import ProfileComparison from './profileComparison';
import PCA from './pca';
import Kataegis from './kataegis';
import Download from './download';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';
import { useSelector, useDispatch } from 'react-redux';
import { actions as visualizationActions } from '../../../services/store/visualization';
import { actions as modalActions } from '../../../services/store/modal';
import {
  value2d,
  filter2d,
  unique2d,
  defaultProfile,
  defaultMatrix,
  defaultFilter,
} from '../../../services/utils';

const actions = { ...visualizationActions, ...modalActions };

export default function Results({ setOpenSidebar }) {
  const dispatch = useDispatch();
  const visualization = useSelector((state) => state.visualization);
  const mergeVisualize = (state) =>
    dispatch(actions.mergeVisualization({ visualize: state }));
  const mergeResults = (state) =>
    dispatch(actions.mergeVisualization({ results: state }));
  const mergeMutationalProfiles = (state) =>
    dispatch(actions.mergeVisualization({ mutationalProfiles: state }));
  const mergeKataegis = (state) =>
    dispatch(actions.mergeVisualization({ kataegis: state }));
  const mergeCosineSimilarity = (state) =>
    dispatch(actions.mergeVisualization({ cosineSimilarity: state }));
  const mergeProfileComparison = (state) =>
    dispatch(actions.mergeVisualization({ profileComparison: state }));
  const mergePCA = (state) =>
    dispatch(actions.mergeVisualization({ pca: state }));
  const mergeError = (msg) =>
    dispatch(actions.mergeModal({ error: { visible: true, message: msg } }));

  const {
    error,
    projectID,
    displayTab,
    svgList,
    matrixList,
  } = visualization.results;

  const { source, loading } = visualization.visualize;
  const mutationalProfiles = visualization.mutationalProfiles;
  const { signatureSetOptions } = visualization.pca;

  // get mapping of plots after retrieving projectID
  useEffect(() => {
    if (source == 'user') {
      if (projectID && !Object.keys(svgList).length) {
        getResults();
      } else if (Object.keys(svgList).length && !signatureSetOptions.length) {
        loadData();
      }
    } else {
      if (
        Object.keys(svgList).length > 0 &&
        !mutationalProfiles.filtered.length
      )
        mapPublicData();
    }
  }, [svgList, projectID]);

  // reload summary information
  async function getResults() {
    const response = await fetch(`api/getResults`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ projectID: projectID }),
    });

    if (response.ok) {
      const {
        svgList,
        statistics,
        matrixList,
        downloads,
      } = await response.json();
      mergeResults({
        svgList: svgList,
        statistics: statistics,
        matrixList: matrixList,
        downloads: downloads,
      });
    } else {
      mergeError(await response.json());
    }
  }

  // retrieve mapping of samples to plots from svgList file
  async function loadData() {
    mergeVisualize({
      loading: {
        active: true,
        content: 'Putting Data Into Session',
        showIndicator: true,
      },
    });

    // Mutational Profiles
    const nameOptions = unique2d('Sample_Name', svgList.columns, svgList.data);
    const selectName = nameOptions[0];
    const filteredPlots = filter2d(selectName, svgList.data);

    const filteredProfileOptions = unique2d(
      'Profile_Type',
      svgList.columns,
      filteredPlots
    );

    const profile = defaultProfile(filteredProfileOptions);

    const filteredMatrixOptions = unique2d(
      'Matrix_Size',
      svgList.columns,
      filter2d(profile, filteredPlots)
    );

    const matrix = defaultMatrix(profile, filteredMatrixOptions);

    const filteredFilterOptions = unique2d(
      'Filter',
      svgList.columns,
      filter2d(matrix, filteredPlots)
    );

    const filter = defaultFilter(filteredFilterOptions);

    const filteredMatrixList = unique2d(
      'Matrix_Size',
      matrixList.columns,
      filter2d(profile, matrixList.data)
    );

    mergeMutationalProfiles({
      filtered: filteredPlots,
      nameOptions: nameOptions,
      profileOptions: filteredProfileOptions,
      matrixOptions: filteredMatrixOptions,
      filterOptions: filteredFilterOptions,
      selectName: selectName,
      selectProfile: profile,
      selectMatrix: matrix,
      selectFilter: filter,
    });

    // Cosine Similarity - Profile Comparison - PCA - Kataegis
    const sampleNameOptions = [
      ...new Set(
        svgList.data.map((row) => {
          if (value2d(row, 'Filter', svgList.columns) != 'NA')
            return `${value2d(row, 'Sample_Name', svgList.columns)}@${value2d(
              row,
              'Filter',
              svgList.columns
            )}`;
          else return value2d(row, 'Sample_Name', svgList.columns);
        })
      ),
    ];
    const profileOptions = unique2d(
      'Profile_Type',
      svgList.columns,
      svgList.data
    );

    const selectProfile = defaultProfile(profileOptions);
    const selectMatrix = defaultMatrix(selectProfile, filteredMatrixOptions);

    const refSignatureSetOptions = await (
      await getRefSigOptions(selectProfile)
    ).json();

    mergeCosineSimilarity({
      withinProfileType: selectProfile,
      refProfileType: selectProfile,
      refSignatureSet: refSignatureSetOptions[0],
      refSignatureSetOptions: refSignatureSetOptions,
      withinMatrixSize: selectMatrix,
      withinMatrixOptions: filteredMatrixList,
      userProfileType: selectProfile,
      userMatrixSize: selectMatrix,
      userMatrixOptions: filteredMatrixOptions,
    });

    mergeProfileComparison({
      withinProfileType: selectProfile,
      withinSampleName1: sampleNameOptions[0],
      withinSampleName2: sampleNameOptions[1],
      sampleOptions: sampleNameOptions,
      refProfileType: selectProfile,
      refSampleName: sampleNameOptions[0],
      refSignatureSet: refSignatureSetOptions[0],
      refSignatureSetOptions: refSignatureSetOptions,
      userProfileType: selectProfile,
      userMatrixSize: selectMatrix,
      userMatrixOptions: filteredMatrixOptions,
      userSampleName: sampleNameOptions[0],
    });

    mergePCA({
      profileType: selectProfile,
      signatureSet: refSignatureSetOptions[0],
      signatureSetOptions: refSignatureSetOptions,
      userProfileType: selectProfile,
      userMatrixSize: selectMatrix,
      userMatrixOptions: filteredMatrixOptions,
    });

    mergeKataegis({
      sample: sampleNameOptions[0],
      sampleOptions: sampleNameOptions,
    });

    mergeVisualize({
      loading: {
        active: false,
      },
    });
    setOpenSidebar(false);
  }

  // retrieve mapping of samples to plots from svgList file
  async function mapPublicData() {
    mergeVisualize({
      loading: {
        active: true,
        content: 'Putting Public Data Into Session',
        showIndicator: true,
      },
    });

    // Mutational Profiles
    const nameOptions = unique2d('Sample', svgList.columns, svgList.data);
    const selectName = mutationalProfiles.selectName || nameOptions[0];
    const profileOptions = [
      ...new Set(
        svgList.data.map(
          (row) => value2d(row, 'Profile', svgList.columns).match(/[a-z]+/gi)[0]
        )
      ),
    ];
    const profile = defaultProfile(profileOptions);
    const matrixOptions = [
      ...new Set(
        svgList.data.map(
          (row) => value2d(row, 'Profile', svgList.columns).match(/\d+/gi)[0]
        )
      ),
    ];
    const matrix = defaultMatrix(profile, matrixOptions);
    const filterOptions = ['NA'];
    const selectFilter = filterOptions[0];

    const filteredPlots = filter2d(
      [selectName, `${profile + matrix}`],
      svgList.data
    );

    const filteredProfileOptions = [
      ...new Set(
        filter2d(selectName, svgList.data).map(
          (row) => value2d(row, 'Profile', svgList.columns).match(/[a-z]+/gi)[0]
        )
      ),
    ].sort();

    const filteredMatrixOptions = [
      ...new Set(
        svgList.data
          .filter(
            (row) =>
              row.includes(selectName) &&
              value2d(row, 'Profile', svgList.columns).indexOf(profile) > -1
          )
          .map(
            (row) => value2d(row, 'Profile', svgList.columns).match(/\d+/gi)[0]
          )
      ),
    ].sort((a, b) => a - b);

    mergeMutationalProfiles({
      filtered: filteredPlots,
      nameOptions: nameOptions,
      profileOptions: filteredProfileOptions,
      matrixOptions: filteredMatrixOptions,
      filterOptions: filterOptions,
      selectName: selectName,
      selectProfile: profile,
      selectMatrix: matrix,
      selectFilter: selectFilter,
    });

    // Cosine Similarity - Profile Comparison - PCA
    const selectProfile = defaultProfile(profileOptions);
    const selectMatrix = defaultMatrix(selectProfile, filteredMatrixOptions);

    const refSignatureSetOptions = await (
      await getRefSigOptions(selectProfile)
    ).json();

    mergeCosineSimilarity({
      withinProfileType: selectProfile,
      refProfileType: selectProfile,
      refSignatureSet: refSignatureSetOptions[0],
      refSignatureSetOptions: refSignatureSetOptions,
      withinMatrixSize: selectMatrix,
      withinMatrixOptions: filteredMatrixOptions,
    });

    mergeProfileComparison({
      withinProfileType: selectProfile,
      withinSampleName1: nameOptions[0],
      withinSampleName2: nameOptions[1],
      sampleOptions: nameOptions,
      refProfileType: selectProfile,
      refSampleName: nameOptions[0],
      refSignatureSet: refSignatureSetOptions[0],
      refSignatureSetOptions: refSignatureSetOptions,
    });

    mergePCA({
      profileType: selectProfile,
      signatureSet: refSignatureSetOptions[0],
      signatureSetOptions: refSignatureSetOptions,
    });

    mergeVisualize({
      loading: {
        active: false,
      },
    });
    setOpenSidebar(false);
  }

  function submitR(fn, args) {
    return fetch(`api/visualizationCalc`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ fn: fn, args: args, projectID: projectID }),
    });
  }

  function getRefSigOptions(profileType) {
    return fetch(`api/visualizationData`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        fn: 'getReferenceSignatureSets',
        args: { profileType: profileType },
      }),
    });
  }

  const examples = [
    { title: 'VCF Example', path: 'vcfExample' },
    { title: 'PCAWG Lung-AdenoCA', path: 'pcawg-lungadenoca' },
    { title: 'PCAWG Lung-SCC', path: 'pcawg-lungscc' },
    { title: 'PCAWG Breast-AdenoCA', path: 'pcawg-breastadenoca' },
    { title: 'PCAWG Skin-Melanoma', path: 'pcawg-skinmelanoma' },
    { title: 'PCAWG PanCancer', path: 'pcawg-pancancer' },
    { title: 'TCGA PanCancer', path: 'tcga-pancancer' },
    {
      title: 'MBD4 defect is associated with hypermutated CpG>TpG pattern',
      external: {
        name: 'PMID: 29760383',
        href: 'https://pubmed.ncbi.nlm.nih.gov/29760383/',
      },
      path: 'mbd4_defect',
    },
  ];

  return (
    <div>
      <LoadingOverlay
        active={loading.active}
        content={loading.content}
        showIndicator={loading.showIndicator}
      />
      {error.length ? (
        <div
          className="border rounded bg-white p-3"
          style={{ minHeight: '420px' }}
        >
          <h4 className="text-danger">{error}</h4>
        </div>
      ) : mutationalProfiles.filtered.length ? (
        <div>
          <div
            style={{
              display: displayTab == 'profilerSummary' ? 'block' : 'none',
              minHeight: '420px',
            }}
          >
            <ProfilerSummary submitR={(fn, args) => submitR(fn, args)} />
          </div>
          <div
            style={{
              display: displayTab == 'mutationalProfiles' ? 'block' : 'none',
              minHeight: '420px',
            }}
          >
            <MutationalProfiles />
          </div>
          <div
            style={{
              display: displayTab == 'cosineSimilarity' ? 'block' : 'none',
              minHeight: '420px',
            }}
          >
            <CosineSimilarity
              getRefSigOptions={(profileType) => getRefSigOptions(profileType)}
              submitR={(fn, args) => submitR(fn, args)}
            />
          </div>
          <div
            style={{
              display: displayTab == 'mutationalPattern' ? 'block' : 'none',
              minHeight: '420px',
            }}
          >
            <MutationalPattern submitR={(fn, args) => submitR(fn, args)} />
          </div>

          <div
            style={{
              display: displayTab == 'profileComparison' ? 'block' : 'none',
              minHeight: '420px',
            }}
          >
            <ProfileComparison
              getRefSigOptions={(profileType) => getRefSigOptions(profileType)}
              submitR={(fn, args) => submitR(fn, args)}
            />
          </div>
          <div
            style={{
              display: displayTab == 'pca' ? 'block' : 'none',
              minHeight: '420px',
            }}
          >
            <PCA
              getRefSigOptions={(profileType) => getRefSigOptions(profileType)}
              submitR={(fn, args) => submitR(fn, args)}
            />
          </div>
          <div
            style={{
              display: displayTab == 'kataegis' ? 'block' : 'none',
              minHeight: '420px',
            }}
          >
            <Kataegis submitR={(fn, args) => submitR(fn, args)} />
          </div>
          <div
            style={{
              display: displayTab == 'download' ? 'block' : 'none',
              minHeight: '420px',
            }}
          >
            <Download />
          </div>
        </div>
      ) : (
        <div
          className="border rounded bg-white p-3"
          style={{ minHeight: '420px' }}
        >
          <h4>Instructions</h4>
          <p>
            Choose a Data Source and its associated options to submit a query
            using the panel on the left
          </p>
          <hr />
          <h4>Data Source</h4>
          <p>Public: Perform analysis using data available on the website</p>
          <p>User: Upload your own data</p>
          <hr />
          <h4>Example Queries</h4>
          {examples.map(({ title, external, path }, index) => (
            <div key={index}>
              <Link to={`/visualization/example/${path}`} disabled>
                <span className="sr-only">{title + ' link'}</span>
                {title}
              </Link>
              {external && (
                <span>
                  {'; '}
                  <a href={external.href} target="_blank">
                    {external.name}
                  </a>
                </span>
              )}
            </div>
          ))}
        </div>
      )}
    </div>
  );
}
