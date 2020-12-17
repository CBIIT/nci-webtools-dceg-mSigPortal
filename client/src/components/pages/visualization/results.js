import React, { useEffect } from 'react';
import { Link } from 'react-router-dom';
import { Card, Nav } from 'react-bootstrap';
import { useSelector } from 'react-redux';
import {
  dispatchError,
  dispatchVisualize,
  dispatchVisualizeResults,
  dispatchMutationalProfiles,
  dispatchCosineSimilarity,
  dispatchProfileComparison,
  dispatchPCA,
} from '../../../services/store';
import ProfilerSummary from './profilerSummary';
import MutationalProfiles from './mutationalProfiles';
import CosineSimilarity from './cosineSimilarity';
import MutationalPattern from './mutationalPattern';
import ProfileComparison from './profileComparison';
import PCA from './pca';
import Download from './download';
import { LoadingOverlay } from '../../controls/loading-overlay/loading-overlay';

const { Header, Body } = Card;
const { Item, Link: CardLink } = Nav;

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

  // retrieve mapping of samples to plots from svgList file
  async function loadData() {
    dispatchVisualize({
      loading: {
        active: true,
        content: 'Putting Data Into Session',
        showIndicator: true,
      },
    });

    const nameOptions = [
      ...new Set(svgList.map(({ Sample_Name }) => Sample_Name)),
    ];
    const profileComparisonSamples = [
      ...new Set(
        svgList.map((plot) => {
          if (plot.Filter != 'NA') return `${plot.Sample_Name}@${plot.Filter}`;
          else return plot.Sample_Name;
        })
      ),
    ];
    const profileOptions = [
      ...new Set(svgList.map((plot) => plot.Profile_Type)),
    ];
    const matrixOptions = [...new Set(svgList.map((plot) => plot.Matrix_Size))];
    const filterOptions = [...new Set(svgList.map((plot) => plot.Filter))];

    const selectName = mutationalProfiles.selectName || nameOptions[0];
    const selectProfile = mutationalProfiles.selectProfile || profileOptions[0];
    const selectMatrix = mutationalProfiles.selectMatrix || matrixOptions[0];
    const selectFilter = mutationalProfiles.selectFilter || filterOptions[0];

    const filteredPlots = svgList.filter(
      (plot) =>
        plot.Sample_Name == selectName &&
        plot.Profile_Type == selectProfile &&
        plot.Matrix_Size == selectMatrix &&
        plot.Filter == selectFilter
    );

    const filteredProfileOptions = [
      ...new Set(
        svgList
          .filter((plot) => plot.Sample_Name == selectName)
          .map((plot) => plot.Profile_Type)
      ),
    ];
    const filteredMatrixOptions = [
      ...new Set(
        svgList
          .filter((plot) => plot.Profile_Type == selectProfile)
          .map((plot) => plot.Matrix_Size)
      ),
    ];
    const filteredFilterOptions = [
      ...new Set(
        svgList
          .filter((plot) => plot.Matrix_Size == selectMatrix)
          .map((plot) => plot.Filter)
      ),
    ];
    const filteredMatrixList = [
      ...new Set(
        matrixList
          .filter((matrix) => matrix.Profile_Type == selectProfile)
          .map((matrix) => matrix.Matrix_Size)
      ),
    ];

    const refSignatureSetOptions = await (
      await getRefSigOptions(profileOptions[0])
    ).json();

    dispatchMutationalProfiles({
      filtered: filteredPlots,
      nameOptions: nameOptions,
      profileOptions: filteredProfileOptions,
      matrixOptions: filteredMatrixOptions,
      filterOptions: filteredFilterOptions,
      selectName: selectName,
      selectProfile: selectProfile,
      selectMatrix: selectMatrix,
      selectFilter: selectFilter,
    });

    dispatchCosineSimilarity({
      withinProfileType: profileOptions[0],
      refProfileType: profileOptions[0],
      refSignatureSet: refSignatureSetOptions[0],
      refSignatureSetOptions: refSignatureSetOptions,
      withinMatrixSize: filteredMatrixList[0],
      withinMatrixOptions: filteredMatrixList,
      userProfileType: profileOptions[0],
      userMatrixSize: filteredMatrixOptions[0],
      userMatrixOptions: filteredMatrixOptions,
    });

    dispatchProfileComparison({
      withinProfileType: profileOptions[0],
      withinSampleName1: profileComparisonSamples[0],
      withinSampleName2: profileComparisonSamples[1],
      sampleOptions: profileComparisonSamples,
      refProfileType: profileOptions[0],
      refSampleName: profileComparisonSamples[0],
      refSignatureSet: refSignatureSetOptions[0],
      refSignatureSetOptions: refSignatureSetOptions,
      userProfileType: profileOptions[0],
      userMatrixSize: filteredMatrixOptions[0],
      userMatrixOptions: filteredMatrixOptions,
      userSampleName: profileComparisonSamples[0],
    });

    dispatchPCA({
      profileType: profileOptions[0],
      signatureSet: refSignatureSetOptions[0],
      signatureSetOptions: refSignatureSetOptions,
      userProfileType: profileOptions[0],
      userMatrixSize: filteredMatrixList[0],
      userMatrixOptions: filteredMatrixOptions,
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

    const nameOptions = [...new Set(svgList.map((plot) => plot.Sample))];
    const profileOptions = [
      ...new Set(svgList.map((plot) => plot.Profile.match(/[a-z]+/gi)[0])),
    ];
    const matrixOptions = [
      ...new Set(svgList.map((plot) => plot.Profile.match(/\d+/gi)[0])),
    ];
    const filterOptions = ['NA'];

    const selectName = mutationalProfiles.selectName || nameOptions[0];
    const selectProfile = mutationalProfiles.selectProfile || profileOptions[0];
    const selectMatrix = mutationalProfiles.selectMatrix || matrixOptions[0];
    const selectFilter = filterOptions[0];

    const filteredPlots = svgList.filter(
      (plot) =>
        plot.Sample == selectName &&
        plot.Profile.indexOf(selectProfile) > -1 &&
        plot.Profile.indexOf(selectMatrix) > -1
    );

    const filteredProfileOptions = [
      ...new Set(
        svgList
          .filter((plot) => plot.Sample == selectName)
          .map((plot) => plot.Profile.match(/[a-z]+/gi)[0])
      ),
    ];

    const filteredMatrixOptions = [
      ...new Set(
        svgList
          .filter(
            (plot) =>
              plot.Sample == selectName &&
              plot.Profile.indexOf(filteredProfileOptions[0]) > -1
          )
          .map((plot) => plot.Profile.match(/\d+/gi)[0])
      ),
    ];

    const refSignatureSetOptions = await (
      await getRefSigOptions(profileOptions[0])
    ).json();

    dispatchMutationalProfiles({
      filtered: filteredPlots,
      nameOptions: nameOptions,
      profileOptions: filteredProfileOptions,
      matrixOptions: filteredMatrixOptions,
      filterOptions: filterOptions,
      selectName: selectName,
      selectProfile: selectProfile,
      selectMatrix: selectMatrix,
      selectFilter: selectFilter,
    });

    dispatchCosineSimilarity({
      withinProfileType: profileOptions[0],
      refProfileType: profileOptions[0],
      refSignatureSet: refSignatureSetOptions[0],
      refSignatureSetOptions: refSignatureSetOptions,
      withinMatrixSize: filteredMatrixOptions[0],
      withinMatrixOptions: filteredMatrixOptions,
    });

    dispatchProfileComparison({
      withinProfileType: profileOptions[0],
      withinSampleName1: nameOptions[0],
      withinSampleName2: nameOptions[1],
      sampleOptions: nameOptions,
      refProfileType: profileOptions[0],
      refSampleName: nameOptions[0],
      refSignatureSet: refSignatureSetOptions[0],
      refSignatureSetOptions: refSignatureSetOptions,
    });

    dispatchPCA({
      profileType: profileOptions[0],
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

  const links = [
    { title: 'Profiler Summary', id: 'profilerSummary' },
    { title: 'Mutational Profiles', id: 'mutationalProfiles' },
    { title: 'Cosine Similarity', id: 'cosineSimilarity' },
    {
      title: 'Mutational Pattern Enrichment Analysis',
      id: 'mutationalPattern',
    },
    { title: 'Profile Comparison', id: 'profileComparison' },
    { title: 'PCA', id: 'pca' },
    { title: 'Download', id: 'download' },
  ];

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
        <Card>
          <Header>
            <Nav variant="pills" defaultActiveKey="#mutationalProfiles">
              {links.map(({ title, id }) => {
                return (
                  <Item key={id}>
                    <CardLink
                      active={displayTab == id}
                      onClick={() =>
                        dispatchVisualizeResults({ displayTab: id })
                      }
                    >
                      {title}
                    </CardLink>
                  </Item>
                );
              })}
            </Nav>
          </Header>
          <Body
            style={{
              display: displayTab == 'profilerSummary' ? 'block' : 'none',
              minHeight: '420px',
            }}
          >
            <ProfilerSummary submitR={(fn, args) => submitR(fn, args)} />
          </Body>
          <Body
            style={{
              display: displayTab == 'mutationalProfiles' ? 'block' : 'none',
              minHeight: '420px',
            }}
          >
            <MutationalProfiles />
          </Body>
          <Body
            style={{
              display: displayTab == 'cosineSimilarity' ? 'block' : 'none',
              minHeight: '420px',
            }}
          >
            <CosineSimilarity
              getRefSigOptions={(profileType) => getRefSigOptions(profileType)}
              submitR={(fn, args) => submitR(fn, args)}
            />
          </Body>
          <Body
            style={{
              display: displayTab == 'mutationalPattern' ? 'block' : 'none',
              minHeight: '420px',
            }}
          >
            <MutationalPattern submitR={(fn, args) => submitR(fn, args)} />
          </Body>

          <Body
            style={{
              display: displayTab == 'profileComparison' ? 'block' : 'none',
              minHeight: '420px',
            }}
          >
            <ProfileComparison
              getRefSigOptions={(profileType) => getRefSigOptions(profileType)}
              submitR={(fn, args) => submitR(fn, args)}
            />
          </Body>
          <Body
            style={{
              display: displayTab == 'pca' ? 'block' : 'none',
              minHeight: '420px',
            }}
          >
            <PCA
              getRefSigOptions={(profileType) => getRefSigOptions(profileType)}
              submitR={(fn, args) => submitR(fn, args)}
            />
          </Body>
          <Body
            style={{
              display: displayTab == 'download' ? 'block' : 'none',
              minHeight: '420px',
            }}
          >
            <Download />
          </Body>
        </Card>
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
