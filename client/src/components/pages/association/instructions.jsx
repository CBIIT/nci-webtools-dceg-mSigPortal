import React from 'react';
import { Table } from 'react-bootstrap';
import { LoadingOverlay } from '@/components/controls/loading-overlay/loading-overlay';

const exposureMeasures = [
  {
    option: 'Signature_exposure',
    meaning:
      'The estimated exposure assigned to each signature in a sample, on the scale stored for that study. Use this for the worked examples.',
  },
  {
    option: 'Signature_exposure_ratio',
    meaning:
      'The signature exposure divided by the sum of signature exposures for that sample in the selected reference set. This is a fraction, not a percentage or an absolute mutation burden.',
  },
  {
    option: 'Signature_exposure_cat',
    meaning: (
      <>
        A two-category exposure classification: <code>Observed</code> and{' '}
        <code>Not_observed</code>. In the reviewed public implementation,
        exposure greater than 100 is classified as Observed; exposure of 100 or
        less is Not_observed. This operational cutoff does not establish
        biological presence or absence.
      </>
    ),
  },
];

const methodPathways = [
  {
    variables: 'Two numeric variables',
    parametric: 'Pearson correlation',
    nonparametric: 'Spearman rank correlation',
  },
  {
    variables: 'Numeric variable and two categories',
    parametric: 'Welch two-sample t-test',
    nonparametric: 'Wilcoxon rank-sum test',
  },
  {
    variables: 'Numeric variable and more than two categories',
    parametric: 'Welch one-way ANOVA',
    nonparametric: 'Kruskal–Wallis test',
  },
  {
    variables: 'Two categorical variables',
    parametric: 'Chi-square pathway',
    nonparametric: 'Chi-square pathway',
  },
];

const outputColumns = [
  {
    id: 'method',
    output: (
      <>
        <code>method</code> and <code>effectsize</code>
      </>
    ),
    interpretation:
      'Identify the statistical test and the scale of the reported effect. Different tests can return different effect measures.',
  },
  {
    id: 'estimate',
    output: <code>estimate</code>,
    interpretation:
      "Read the effect's magnitude and direction on its stated scale. In a linear model, this is a regression coefficient, not an odds ratio.",
  },
  {
    id: 'n.obs',
    output: <code>n.obs</code>,
    interpretation:
      'The number of analyzed observations where provided. It can be smaller than the study size and can vary by signature. The current multivariable table does not display model N.',
  },
  {
    id: 'p.value',
    output: <code>p.value</code>,
    interpretation: (
      <>
        The unadjusted p-value where returned. A value displayed as{' '}
        <code>0</code> may reflect rounding; do not report it as exactly zero.
      </>
    ),
  },
  {
    id: 'fdr',
    output: (
      <>
        <code>fdr</code> and <code>fdr.method</code>
      </>
    ),
    interpretation:
      'Where supplied, the adjusted p-value and correction method. For the common frequentist univariable tests, BH adjusts across valid signature tests returned for that run. It does not cover additional analyses performed in separate runs.',
  },
  {
    id: 'term',
    output: (
      <>
        <code>term</code> and <code>std.error</code>
      </>
    ),
    interpretation:
      'In multivariable results, the model coefficient label and its standard error. The standard error is not itself a confidence interval.',
  },
];

export default function Instructions({ loading }) {
  return (
    <div className="bg-white border rounded py-3 px-4 guidance">
      <LoadingOverlay active={loading} />
      <h4>Instructions</h4>
      <p>
        The Association module explores relationships between mutational
        signature exposure and clinical or genomic variables in public studies
        available in mSigPortal. <b>Univariable</b> analysis examines one
        association variable at a time. <b>Multivariable</b> analysis models
        signature exposure using several association variables together,
        allowing each association to be evaluated conditional on the other
        selected variables.
      </p>
      <p>
        These are exploratory associations. They do not establish that a
        variable causes a mutational process, or that a mutational process
        causes a clinical or genomic feature.
      </p>
      <hr />
      <h4>Data Source</h4>
      <p>Public: Perform analysis using data available on the website</p>
      {/* <p>User: Upload your own data</p> */}

      <hr />
      <h5>1. Load a study</h5>
      <ol>
        <li>Open the Association module.</li>
        <li>
          In the input panel on the left, select <b>Study</b>,{' '}
          <b>Experimental Strategy</b>, <b>Reference Signature Set</b>, and{' '}
          <b>Cancer Type or Group</b>.
        </li>
        <li>
          Click <b>Load Study</b>. When loading is complete, the{' '}
          <b>Univariable</b> tab opens and the <b>Variables</b> table lists the
          variables available for the selected study.
        </li>
        <li>
          Browse or search the table before choosing variables. The{' '}
          <code>variable_value_type</code> column describes their stored types.
          Available choices depend on the selected study and cohort.
        </li>
      </ol>
      <p>
        The current interface uses public data already available in mSigPortal;
        no user-data upload is needed for these examples. The left input panel
        may collapse after loading. Use <b>Show or hide input sidebar</b> to
        reopen it.
      </p>

      <h5>2. Choose the signature measure</h5>
      <p>
        Select <b>Signature Exposure Variable</b> in the Variables section
        before clicking <b>Load Data</b>.
      </p>
      <Table striped bordered size="sm" responsive className="guidance-table">
        <thead>
          <tr>
            <th>Option</th>
            <th>Meaning</th>
          </tr>
        </thead>
        <tbody>
          {exposureMeasures.map(({ option, meaning }) => (
            <tr key={option}>
              <td>
                <code>{option}</code>
              </td>
              <td>{meaning}</td>
            </tr>
          ))}
        </tbody>
      </Table>
      <p>
        Use a numeric signature measure for the multivariable linear-model
        example. Do not interpret a categorical exposure with <code>lm</code> or
        an unspecified <code>glm</code> as logistic regression.
      </p>

      <h5>3. Run a univariable analysis</h5>
      <ol>
        <li>
          Select the <b>Univariable</b> tab.
        </li>
        <li>
          Choose <b>Variable Source</b>, <b>Data Type</b>, and{' '}
          <b>Variable Name</b>. For example, a gene mutation variable may appear
          under <code>genomic_data</code> and <code>driver_mutation</code>.
        </li>
        <li>
          Set any <b>Variable Filtering</b> options needed for the selected
          association variable. For an introductory analysis, leave{' '}
          <b>Threshold</b> blank, <b>log₂</b> unchecked, and <b>Collapse</b>{' '}
          unchanged.
        </li>
        <li>
          Choose <b>Signature Exposure Variable</b>, then click <b>Load Data</b>
          .
        </li>
        <li>
          In <b>Parameters</b>, review <b>Signature Exposure Filtering</b> and
          select <b>Method</b>. Filtering in this section applies to the
          signature measure; it is separate from filtering of the association
          variable above.
        </li>
        <li>
          Click <b>Calculate</b> and wait for the results table and plot.
        </li>
        <li>
          Review the table, which contains results across signatures. Select{' '}
          <b>Signature Name</b> to display the plot for a particular signature.
          Changing this selector updates the selected-signature output; it does
          not mean the all-signature table is restricted to that signature.
        </li>
      </ol>

      <h6>Filtering and transformations</h6>
      <p>
        For univariable inputs treated as numeric by the analysis,{' '}
        <b>Threshold</b> retains values strictly greater than the entered value.
        Leaving it blank differs from entering zero: a threshold of zero
        excludes zero values. The analysis may treat a stored numeric variable
        with few distinct values as categorical; check the returned method.
      </p>
      <p>
        <b>log₂</b> applies <code>log₂(x)</code> to positive values and excludes
        values of zero or less; it does not apply <code>log₂(x + 1)</code>. The
        number of analyzed samples can therefore differ between signatures.
        Checking log₂ is an analytical choice, not a required setup step. For a
        rank-based test, the transformation preserves the ranking of retained
        positive values, but excluding zero values changes the population being
        analyzed.
      </p>
      <p>
        <b>Collapse</b> groups categories of an association variable. Use it
        only when the resulting grouping matches the scientific question, and
        record the grouping used. Leave it unchanged in the basic examples.
      </p>

      <h6>Choosing a method</h6>
      <p>
        The actual test depends on the variable types and the number of
        categories. The following are common analysis pathways; check the
        returned <code>method</code> field for the test actually run.
      </p>
      <Table striped bordered size="sm" responsive className="guidance-table">
        <thead>
          <tr>
            <th>Variables being compared</th>
            <th>parametric</th>
            <th>nonparametric</th>
          </tr>
        </thead>
        <tbody>
          {methodPathways.map(({ variables, parametric, nonparametric }) => (
            <tr key={variables}>
              <td>{variables}</td>
              <td>{parametric}</td>
              <td>{nonparametric}</td>
            </tr>
          ))}
        </tbody>
      </Table>
      <p>
        Choose a method appropriate to the scientific question, distribution,
        category sizes, and independence of samples. A nonparametric option does
        not resolve sparse categories, confounding, or dependence between
        samples. An overall test across more than two categories does not
        identify which individual pairs differ. Advanced methods require their
        own assumptions and interpretation.
      </p>

      <h5>4. Run a multivariable analysis</h5>
      <ol>
        <li>
          Select the <b>Multivariable</b> tab.
        </li>
        <li>
          Choose <b>Variable Source</b>, <b>Data Type</b>, and{' '}
          <b>Variable Name</b> for the first association variable.
        </li>
        <li>
          Click <b>Add Association Variable</b> to add each additional predictor
          or covariate. Use <b>Remove</b> to remove an unwanted row.
        </li>
        <li>
          Add only the covariates needed for the planned analysis. Missing data,
          sparse categories, and strongly correlated predictors can reduce
          precision or make a model unstable. A categorical variable can
          contribute more than one model coefficient.
        </li>
        <li>
          Choose a numeric <b>Signature Exposure Variable</b>, then click{' '}
          <b>Load Data</b>.
        </li>
        <li>
          For the basic example, leave all thresholds blank, leave all log₂
          boxes unchecked, leave Collapse unchanged, select{' '}
          <b>
            Method = <code>lm</code>
          </b>
          , and click <b>Calculate</b>.
        </li>
        <li>
          Review the model terms and estimates in the results table. Select{' '}
          <b>Signature Name</b> to inspect a signature-specific coefficient
          plot, and verify that the plotted estimates match that signature's
          table rows.
        </li>
      </ol>
      <p>
        The linear model treats the signature measure as the outcome and the
        selected association variables as predictors. Each coefficient describes
        an association conditional on the other included predictors. For a
        numeric predictor, it is the expected change in the outcome per one-unit
        increase in that predictor. For a categorical predictor, it compares a
        coded level with its reference level. Identify the category coding
        before assigning a biological interpretation to a term such as{' '}
        <code>EGFROther</code>.
      </p>
      <p>
        For this multivariable example, leave exposure <b>Threshold</b> blank
        and <b>log₂</b> unchecked. The <code>glm</code> option uses
        Gaussian/identity defaults; do not interpret it as automatic logistic or
        count regression.
      </p>

      <h5>5. Read and save the results</h5>
      <Table striped bordered size="sm" responsive className="guidance-table">
        <thead>
          <tr>
            <th>Output</th>
            <th>How to interpret it</th>
          </tr>
        </thead>
        <tbody>
          {outputColumns.map(({ id, output, interpretation }) => (
            <tr key={id}>
              <td>{output}</td>
              <td>{interpretation}</td>
            </tr>
          ))}
        </tbody>
      </Table>
      <p>
        Covariate adjustment and multiple-testing correction are different.
        Multivariable output contains unadjusted term-level p-values; do not
        describe these as FDR-adjusted. Assess effect sizes and uncertainty
        alongside p-values. A missing result row can indicate insufficient or
        non-estimable data, not evidence of no association.
      </p>
      <p>
        Use <b>Download Association Result</b> for the results table and{' '}
        <b>Download Plot</b> for the selected plot. <b>Download Data</b>{' '}
        provides data associated with the selected signature, but should not be
        assumed to be the final filtered/transformed analysis matrix.{' '}
        <b>Download Variable Data</b> is available above the study variable
        table.
      </p>
      <p>
        Record the study, strategy, signature set, cohort, variables, exposure
        measure, filters, transformations, method, selected signature, category
        coding, and access date. Record the number analyzed when available; seek
        the model N before formal reporting if it is absent from the output.
        Preserve full numerical precision in exported results when available.
      </p>

      <h5>6. Change or troubleshoot an analysis</h5>
      <ul>
        <li>
          <b>Different variables or parameters, same study:</b> Use <b>Reset</b>{' '}
          in the Variables section, recheck every setting, then repeat{' '}
          <b>Load Data</b> and <b>Calculate</b>. Some previous choices may
          remain selected; Reset should not be assumed to restore every default.
        </li>
        <li>
          <b>Different study or cohort:</b> Reopen the left input panel, use its{' '}
          <b>Reset</b>, select the new study settings, and click{' '}
          <b>Load Study</b>.
        </li>
        <li>
          <b>Analysis tabs are disabled:</b> Complete <b>Load Study</b> first.
        </li>
        <li>
          <b>Calculate is not available:</b> Select an association variable and
          a signature exposure variable, click <b>Load Data</b>, and wait for
          loading to finish.
        </li>
        <li>
          <b>Few or no results:</b> Check variable availability, missingness,
          retained sample counts, zero-value exclusions, constant values, and
          small categories. Do not relax filters solely to obtain a small
          p-value.
        </li>
        <li>
          <b>Plot and table disagree:</b> Confirm the selected signature and
          settings. If the discrepancy persists, report the selected settings
          and a screenshot using <b>Contact Us</b>.
        </li>
      </ul>
    </div>
  );
}
