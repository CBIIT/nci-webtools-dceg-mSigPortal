import React, { useEffect } from 'react';
import { Link } from 'react-router-dom';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchVisualize,
  dispatchVisualizeResults,
  dispatchMutationalProfiles,
  dispatchCosineSimilarity,
  dispatchProfileComparison,
  dispatchRainfall,
  dispatchPCA,
} from '../../../services/store';
import ProfilerSummary from './profilerSummary';
import MutationalProfiles from './mutationalProfiles';
import CosineSimilarity from './cosineSimilarity';
import MutationalPattern from './mutationalPattern';
import ProfileComparison from './profileComparison';
import PCA from './pca';
import Rainfall from './rainfall';
import Download from './download';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

export default function Results({ setOpenSidebar }) {
  const { error, projectID, displayTab, svgList, matrixList } = useSelector(
    (state) => state.visualizeResults
  );
  const { source, loading } = useSelector((state) => state.visualize);
  const mutationalProfiles = useSelector((state) => state.mutationalProfiles);
  const { signatureSetOptions } = useSelector((state) => state.pca);

  // get mapping of plots after retrieving projectID
  useEffect(() => {
    if (source == 'user') {
      if (projectID && !svgList.length) {
        getResultData();
      } else if (svgList.length && !signatureSetOptions.length) loadData();
    } else {
      if (svgList.length > 0 && !mutationalProfiles.filtered.length)
        mapPublicData();
    }
  }, [svgList, projectID, source]);

  // reload summary information
  async function getResultData() {
    const response = await fetch(`api/getResultData`, {
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
      dispatchVisualizeResults({
        svgList: svgList,
        statistics: statistics,
        matrixList: matrixList,
        downloads: downloads,
      });
    } else {
      dispatchError(await response.json());
    }
  }

  function defaultProfile(profileOptions) {
    if (profileOptions.includes('SBS')) return 'SBS';
    if (profileOptions.includes('DBS')) return 'DBS';
    if (profileOptions.includes('ID')) return 'ID';
  }

  function defaultMatrix(profile, matrixOptions) {
    if (profile == 'SBS')
      return matrixOptions.includes('96') ? '96' : matrixOptions[0];

    if (profile == 'DBS')
      return matrixOptions.includes('78') ? '78' : matrixOptions[0];

    if (profile == 'ID')
      return matrixOptions.includes('83') ? '83' : matrixOptions[0];
  }

  function defaultFilter(filterOptions) {
    return filterOptions.includes('NA') ? 'NA' : filterOptions[0];
  }

  // retrieve mapping of samples to plots from svgList file
  async function loadData() {
    dispatchVisualize({
      loading: {
        active: true,
        content: 'Putting Data Into Session',
        showIndicator: true,
      },
    });

    // Mutational Profiles
    const nameOptions = [
      ...new Set(svgList.map(({ Sample_Name }) => Sample_Name)),
    ];
    const selectName = nameOptions[0];

    const filteredPlots = svgList.filter(
      (plot) => plot.Sample_Name == selectName
    );

    const filteredProfileOptions = [
      ...new Set(filteredPlots.map((plot) => plot.Profile_Type)),
    ].sort();
    const mpProfile = defaultProfile(filteredProfileOptions);

    const filteredMatrixOptions = [
      ...new Set(
        filteredPlots
          .filter((plot) => plot.Profile_Type == mpProfile)
          .map((plot) => plot.Matrix_Size)
      ),
    ].sort((a, b) => a - b);
    const mpMatrix = defaultMatrix(mpProfile, filteredMatrixOptions);

    const filteredFilterOptions = [
      ...new Set(
        filteredPlots
          .filter((plot) => plot.Matrix_Size == mpMatrix)
          .map((plot) => plot.Filter)
      ),
    ];
    const mpFilter = defaultFilter(filteredFilterOptions);

    const filteredMatrixList = [
      ...new Set(
        matrixList
          .filter((matrix) => matrix.Profile_Type == mpProfile)
          .map((matrix) => matrix.Matrix_Size)
      ),
    ].sort((a, b) => a - b);

    dispatchMutationalProfiles({
      filtered: filteredPlots,
      nameOptions: nameOptions,
      profileOptions: filteredProfileOptions,
      matrixOptions: filteredMatrixOptions,
      filterOptions: filteredFilterOptions,
      selectName: selectName,
      selectProfile: mpProfile,
      selectMatrix: mpMatrix,
      selectFilter: mpFilter,
    });

    // Cosine Similarity - Profile Comparison - PCA - Rainfall
    const sampleNameOptions = [
      ...new Set(
        svgList.map((plot) => {
          if (plot.Filter != 'NA') return `${plot.Sample_Name}@${plot.Filter}`;
          else return plot.Sample_Name;
        })
      ),
    ];
    const profileOptions = [
      ...new Set(svgList.map((plot) => plot.Profile_Type)),
    ].sort();

    const selectProfile = defaultProfile(profileOptions);
    const selectMatrix = defaultMatrix(selectProfile, filteredMatrixOptions);

    const refSignatureSetOptions = await (
      await getRefSigOptions(selectProfile)
    ).json();

    dispatchCosineSimilarity({
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

    dispatchProfileComparison({
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

    dispatchPCA({
      profileType: selectProfile,
      signatureSet: refSignatureSetOptions[0],
      signatureSetOptions: refSignatureSetOptions,
      userProfileType: selectProfile,
      userMatrixSize: selectMatrix,
      userMatrixOptions: filteredMatrixOptions,
    });

    dispatchRainfall({
      sample: sampleNameOptions[0],
      sampleOptions: sampleNameOptions,
    });

    dispatchVisualize({
      loading: {
        active: false,
      },
    });
    setOpenSidebar(false);
  }

  // retrieve mapping of samples to plots from svgList file
  async function mapPublicData() {
    dispatchVisualize({
      loading: {
        active: true,
        content: 'Putting Public Data Into Session',
        showIndicator: true,
      },
    });

    // Mutational Profiles
    const nameOptions = [...new Set(svgList.map((plot) => plot.Sample))];
    const selectName = mutationalProfiles.selectName || nameOptions[0];

    const profileOptions = [
      ...new Set(svgList.map((plot) => plot.Profile.match(/[a-z]+/gi)[0])),
    ];

    const matrixOptions = [
      ...new Set(svgList.map((plot) => plot.Profile.match(/\d+/gi)[0])),
    ];

    const filterOptions = ['NA'];
    const selectFilter = filterOptions[0];

    const filteredPlots = svgList.filter(
      (plot) =>
        plot.Sample == selectName &&
        plot.Profile.indexOf(defaultProfile(profileOptions)) > -1 &&
        plot.Profile.indexOf(
          defaultMatrix(defaultProfile(profileOptions), matrixOptions)
        ) > -1
    );

    const filteredProfileOptions = [
      ...new Set(
        svgList
          .filter((plot) => plot.Sample == selectName)
          .map((plot) => plot.Profile.match(/[a-z]+/gi)[0])
      ),
    ].sort();
    const mpProfile =
      mutationalProfiles.selectProfile ||
      defaultProfile(filteredProfileOptions);

    const filteredMatrixOptions = [
      ...new Set(
        svgList
          .filter(
            (plot) =>
              plot.Sample == selectName && plot.Profile.indexOf(mpProfile) > -1
          )
          .map((plot) => plot.Profile.match(/\d+/gi)[0])
      ),
    ].sort((a, b) => a - b);
    const mpMatrix =
      mutationalProfiles.selectMatrix ||
      defaultMatrix(mpProfile, filteredMatrixOptions);

    dispatchMutationalProfiles({
      filtered: filteredPlots,
      nameOptions: nameOptions,
      profileOptions: filteredProfileOptions,
      matrixOptions: filteredMatrixOptions,
      filterOptions: filterOptions,
      selectName: selectName,
      selectProfile: mpProfile,
      selectMatrix: mpMatrix,
      selectFilter: selectFilter,
    });

    // Cosine Similarity - Profile Comparison - PCA
    const selectProfile = defaultProfile(profileOptions);
    const selectMatrix = defaultMatrix(selectProfile, filteredMatrixOptions);

    const refSignatureSetOptions = await (
      await getRefSigOptions(selectProfile)
    ).json();

    dispatchCosineSimilarity({
      withinProfileType: selectProfile,
      refProfileType: selectProfile,
      refSignatureSet: refSignatureSetOptions[0],
      refSignatureSetOptions: refSignatureSetOptions,
      withinMatrixSize: selectMatrix,
      withinMatrixOptions: filteredMatrixOptions,
    });

    dispatchProfileComparison({
      withinProfileType: selectProfile,
      withinSampleName1: nameOptions[0],
      withinSampleName2: nameOptions[1],
      sampleOptions: nameOptions,
      refProfileType: selectProfile,
      refSampleName: nameOptions[0],
      refSignatureSet: refSignatureSetOptions[0],
      refSignatureSetOptions: refSignatureSetOptions,
    });

    dispatchPCA({
      profileType: selectProfile,
      signatureSet: refSignatureSetOptions[0],
      signatureSetOptions: refSignatureSetOptions,
    });

    dispatchVisualize({
      loading: {
        active: false,
      },
    });
    setOpenSidebar(false);
  }

  function submitR(fn, args) {
    return fetch(`api/visualizeR`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ fn: fn, args: args, projectID: projectID }),
    });
  }

  function getRefSigOptions(profileType) {
    return fetch(`api/getReferenceSignatureSets`, {
      method: 'POST',
      headers: {
        Accept: 'application/json',
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ profileType: profileType }),
    });
  }

  const examples = [
    { title: 'VCF Example', path: 'vcfExample' },
    {
      title: 'MBD4 defect is associated with hypermutated CpG>TpG pattern',
      external: {
        name: 'PMID: 29760383',
        href: 'https://pubmed.ncbi.nlm.nih.gov/29760383/',
      },
      path: 'tcgaPanCancer',
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
            <MutationalProfiles
              defaultMatrix={(profile, matrixOptions) =>
                defaultMatrix(profile, matrixOptions)
              }
            />
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
              defaultMatrix={(profile, matrixOptions) =>
                defaultMatrix(profile, matrixOptions)
              }
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
              defaultMatrix={(profile, matrixOptions) =>
                defaultMatrix(profile, matrixOptions)
              }
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
              defaultMatrix={(profile, matrixOptions) =>
                defaultMatrix(profile, matrixOptions)
              }
            />
          </div>
          <div
            style={{
              display: displayTab == 'rainfall' ? 'block' : 'none',
              minHeight: '420px',
            }}
          >
            <Rainfall submitR={(fn, args) => submitR(fn, args)} />
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
          <h2>Instructions</h2>
          <hr />
          <h4>Data Souce</h4>
          <p>Public: Perform analysis using data available on the website</p>
          <p>User: Upload Sample Variants and specify parameters </p>
          <hr />
          <h4>Examples Queries</h4>
          {examples.map(({ title, external, path }, index) => (
            <div key={index}>
              <Link to={`visualization/example/${path}`}>
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
