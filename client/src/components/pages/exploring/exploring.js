import React from 'react';
import { Card, Nav } from 'react-bootstrap';
// import ReferenceSignatures from './referenceSignatures';

const { Header, Body } = Card;
const { Item, Link } = Nav;

export default function Explore() {
  return (
    <div className="position-relative">
      {/* <Card>
        <Header>
          <Nav variant="pills" defaultActiveKey="#mutationalProfiles">
            {[
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
            ].map(({ title, id }) => {
              return (
                <Item key={id}>
                  <Link
                    active={displayTab == id}
                    onClick={() => dispatchVisualizeResults({ displayTab: id })}
                  >
                    {title}
                  </Link>
                </Item>
              );
            })}
          </Nav>
        </Header>
        <Body
          style={{
            display: displayTab == 'profilerSummary' ? 'block' : 'none',
          }}
        >
          <ProfilerSummary submitR={(fn, args) => submitR(fn, args)} />
        </Body>
        <Body
          style={{
            display: displayTab == 'mutationalProfiles' ? 'block' : 'none',
          }}
        >
          <MutationalProfiles />
        </Body>
        <Body
          style={{
            display: displayTab == 'cosineSimilarity' ? 'block' : 'none',
          }}
        >
          <CosineSimilarity
            getRefSigOptions={(profileType) => getRefSigOptions(profileType)}
            downloadResults={(path) => downloadResults(path)}
            submitR={(fn, args) => submitR(fn, args)}
          />
        </Body>
        <Body
          style={{
            display: displayTab == 'mutationalPattern' ? 'block' : 'none',
          }}
        >
          <MutationalPattern
            downloadResults={(path) => downloadResults(path)}
            submitR={(fn, args) => submitR(fn, args)}
          />
        </Body>

        <Body
          style={{
            display: displayTab == 'profileComparison' ? 'block' : 'none',
          }}
        >
          <ProfileComparison
            getRefSigOptions={(profileType) => getRefSigOptions(profileType)}
            submitR={(fn, args) => submitR(fn, args)}
          />
        </Body>
        <Body style={{ display: displayTab == 'pca' ? 'block' : 'none' }}>
          <PCA
            getRefSigOptions={(profileType) => getRefSigOptions(profileType)}
            downloadResults={(path) => downloadResults(path)}
            submitR={(fn, args) => submitR(fn, args)}
          />
        </Body>
        <Body style={{ display: displayTab == 'download' ? 'block' : 'none' }}>
          <Download />
        </Body>
      </Card> */}
    </div>
  );
}
