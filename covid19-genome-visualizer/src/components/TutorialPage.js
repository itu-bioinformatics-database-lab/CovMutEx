import React from 'react';

const TutorialPage = () => {
  // A helper object for text styles to reduce repetition, though using them inline is also perfectly fine.
  const sectionHeading = "text-2xl font-bold text-slate-800 mt-10 mb-4 pb-2 border-b border-slate-200";
  const subHeading = "text-xl font-semibold text-slate-700 mt-8 mb-3";
  const bodyText = "mb-5";
  const listStyles = "list-disc list-outside pl-5 mb-5 space-y-3";
  const strongHighlight = "font-semibold text-blue-600";

  // A small component for the color swatch to keep the JSX cleaner
  const ColorSwatch = ({ color, name }) => (
    <li className="flex items-center">
      <span className={`inline-block w-3.5 h-3.5 mr-2 border border-slate-300 rounded-sm ${color}`}></span>
      {name}
    </li>
  );

  return (
    <article className="max-w-6xl mx-auto my-8 p-6 sm:p-8 bg-white rounded-lg shadow-md font-sans text-slate-700 leading-relaxed">
      <h1 className="text-4xl font-bold text-slate-900 mb-6 pb-3 border-b-2 border-slate-200">
        Covid Mutation Explorer Tutorial
      </h1>

      <p className={bodyText}>
        This tutorial introduces Covid Mutation Explorer, a powerful, web-based platform designed to predict and visualize future mutations of SARS-CoV-2. Covid Mutation Explorer equips researchers studying viral evolution and users curious about how the virus changes over time with intuitive tools to explore and anticipate mutation patterns. This guide will walk through its functionality and unique features.
      </p>

      <h2 className={sectionHeading}>Getting Started: Your First Prediction</h2>

      <p className={bodyText}>
        To begin, select a prediction model. Covid Mutation Explorer offers three powerful options, each tailored to different research needs:
      </p>
      <ul className={listStyles}>
        <li><strong>The Single Input Ensemble Model:</strong> Combines several prediction algorithms using one primary data source.</li>
        <li><strong>The Balanced Model:</strong> Designed to reduce bias by equally representing various viral variants.</li>
        <li><strong>The Multi-Input Ensemble Model:</strong> Draws from multiple data sources and prediction techniques for enhanced accuracy.</li>
      </ul>
      <p className={bodyText}>
        This demonstration uses the <strong className={strongHighlight}>Balanced Model</strong>.
      </p>

      <p className={bodyText}>
        Next, choose a variant ID, which corresponds to a specific SARS-CoV-2 variant sequenced in the GISAID database. Then, input the number of elapsed days to tell Covid Mutation Explorer how far into the future to project potential mutations from the date the variant first appeared.
      </p>
      <p className={bodyText}>
        A dropdown menu for protein regions is also available, representing parts of the viral genome like the spike protein. For the first run, this will remain unselected to analyze predictions across the entire genome. Click <strong className={strongHighlight}>Predict</strong>.
      </p>

      <h2 className={sectionHeading}>Understanding the Results</h2>
      <p className={bodyText}>
        Once the prediction runs, results are organized across three panels, each offering a different layer of insight:
      </p>
      <ul className={listStyles}>
        <li>
          <strong>On the left</strong>, the protein region interface allows you to toggle annotations and quickly navigate different parts of the genome. Hovering over a region highlights it in the main view.
        </li>
        <li>
          <strong>On the right</strong>, an interactive pie chart visualizes the probabilities of mutations. Hover to get detailed info, or click <strong className={strongHighlight}>Normalize</strong> to adjust for coding region length differences.
        </li>
        <li>
          <strong>At the center</strong> is the main visualization panel—a display of mutation probabilities as stacked bar charts. Each nucleotide is color-coded:
          <ul className="list-none mt-3 space-y-2">
            <ColorSwatch color="bg-red-500" name="Red for A" />
            <ColorSwatch color="bg-green-500" name="Green for T" />
            <ColorSwatch color="bg-yellow-400" name="Yellow for G" />
            <ColorSwatch color="bg-blue-500" name="Blue for C" />
          </ul>
          The Y-axis uses a logarithmic scale, making it easy to compare both common and rare mutations side by side.
        </li>
      </ul>

      <h2 className={sectionHeading}>Interacting with the Visualization</h2>
      
      <h3 className={subHeading}>Multi-Level Zoom</h3>
      <p className={bodyText}>
        Covid Mutation Explorer features an intuitive, multi-level zoom system. Use the plus and minus buttons—or your mouse wheel—to zoom in for a more granular view. Initially, each bar represents 25 nucleotide positions. Zooming in further reveals an intermediate, high-resolution view to analyze the genome at the nucleotide level.
      </p>
      
      <h3 className={subHeading}>WebLogo Mode</h3>
      <p className={bodyText}>
        When fully zoomed in, the platform automatically switches to WebLogo mode—a visual sequence logo that shows conservation at each position. Here, taller letters mean high conservation, while mixed letters reveal variation. The reference sequence is also visible, giving a clear sense of how each position is evolving.
      </p>

      <h3 className={subHeading}>Navigation and Tooltips</h3>
      <p className={bodyText}>
        To navigate the genome, use the pan controls at the bottom or click and drag left or right. Even in WebLogo mode, panning occurs seamlessly without losing detail, which is ideal for exploring mutations across neighboring regions.
      </p>
      <p className={bodyText}>
        Hovering over any bar in the chart brings up a detailed tooltip, showing:
      </p>
      <ul className={listStyles}>
        <li>The exact genomic position</li>
        <li>The reference nucleotide</li>
        <li>Probabilities for each possible mutation</li>
        <li>The total mutation probability</li>
      </ul>
      <p className={bodyText}>
        The stacked bar format makes it easy to compare the likelihood of different substitutions at a glance.
      </p>

      <h2 className={sectionHeading}>Example: Region-Specific Prediction</h2>
      <p className={bodyText}>
        The following example demonstrates a region-specific prediction. Using the <strong className={strongHighlight}>Balanced Model</strong> again, select a specific protein region—<strong className={strongHighlight}>ORF1AB</strong>. After clicking <strong className={strongHighlight}>Predict</strong>, the pie chart adapts to show data only within the selected region. In the main panel, mutation predictions are limited to <strong className={strongHighlight}>ORF1AB</strong>, and from here, it is possible to zoom in to explore the WebLogo breakdown for this specific protein region.
      </p>

      <h2 className={sectionHeading}>Learn More</h2>
      <p className={bodyText}>
        This guide has covered the core functionality of Covid Mutation Explorer, from setting up predictions to exploring rich, interactive visualizations.
      </p>
      <p className={bodyText}>
        For more information on the methodology behind the prediction models or the data sources used, refer to the <strong className={strongHighlight}>About</strong> page for an in-depth explanation. To inquire about how Covid Mutation Explorer can support specific research, use the <strong className={strongHighlight}>Contact</strong> page. Covid Mutation Explorer is a valuable ally in the exploration of viral mutation dynamics.
      </p>
    </article>
  );
};

export default TutorialPage;
