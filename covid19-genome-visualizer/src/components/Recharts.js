import React, { useEffect, useMemo, useRef, useState } from "react";
import { useSelector } from "react-redux";
import Chart from "chart.js/auto";
import annotationPlugin from "chartjs-plugin-annotation";
import zoomPlugin from "chartjs-plugin-zoom";
import { getColorForNucleotide, nucleotides } from "../helpers/helperFunctions";
import { proteinRegions } from "../data/proteinRegions";
import {
  proteinRegionColorMap,
  proteinRegionColorMapAnnotations,
} from "../utils/proteinRegionColorMap";
import SidePanel from "./SidePanel";
import Nav from "./Nav";

Chart.register(annotationPlugin, zoomPlugin);

let chunkView = true;

const GenomeChart = ({ genomeData, genomeSequence }) => {
  const chartRef = useRef(null);

  const [activeProtein, setActiveProtein] = useState(null);
  const [showFullAnnotation, setShowFullAnnotation] = useState(false);
  const [hoveredPosition, setHoveredPosition] = useState(null);

  const selectedProteinRegion = useSelector(
    (state) => state.genome.selectedProteinRegion
  );
  const [zoomLevel, setZoomLevel] = useState(selectedProteinRegion ? 1 : 25);
  const genomeState = useSelector((state) => state.genome);
  const zoomThreshold = 1000;
  let [decimateFactor, setDecimateFcator] = useState(
    selectedProteinRegion || zoomLevel >= zoomThreshold ? 1 : 25
  );

  decimateFactor = 25;
  useEffect(() => {
    if (selectedProteinRegion) {
      setDecimateFcator(1); // Use full-resolution data for protein regions
      console.log("PR SELECTED DECIMATED DATA", decimateFactor);
    } else if (zoomLevel >= zoomThreshold) {
      setDecimateFcator(1); // Full-resolution data for high zoom levels
    } else {
      setDecimateFcator(25); // Decimated data for lower zoom levels
    }
  }, [selectedProteinRegion, zoomLevel]);

  // console.log(selectedProteinRegion, "ADFFFFFFFFFFFFs");
  // console.log(zoomLevel >= zoomThreshold);
  // console.log(zoomLevel, "ZOOOMMMM");
  // console.log("Decimate Factor", decimateFactor);
  // console.log(genomeSequence[29042], "29028 position");
  let currentDecimateFactor =
    selectedProteinRegion && proteinRegions[selectedProteinRegion]
      ? selectedProteinRegion === "ORF1ab"
        ? zoomLevel >= zoomThreshold
          ? 1
          : 25
        : 1
      : zoomLevel >= zoomThreshold
      ? 1
      : decimateFactor;

  const decimatedData = useMemo(() => {
    if (!genomeData || genomeData.length === 0) return [];

    return nucleotides.map((nucleotide, nucleotideIndex) => {
      return genomeData[nucleotideIndex].reduce((acc, curr, index) => {
        if (index % decimateFactor === 0) {
          const slice = genomeData[nucleotideIndex].slice(
            index,
            index + decimateFactor
          );
          const sum = slice.reduce((sum, value) => sum + value, 0);
          acc.push(sum / decimateFactor);
        }
        return acc;
      }, []);
    });
  }, [genomeData, decimateFactor]);

  const getMaxValue = (data) => {
    if (!data || data.length === 0) return 1;
    return Math.max(...data.flat());
  };

  const maxValue = useMemo(() => getMaxValue(decimatedData), [decimatedData]);

  const createAnnotations = () => {
    if (activeProtein) {
      const range = proteinRegions[activeProtein];
      const [start, end] = range
        .split("-")
        .map((e) => Math.floor(e / decimateFactor));
      return [
        {
          display: true,
          type: "box",
          xMin: start,
          xMax: end,
          yMin: 0,
          yMax: 3,
          backgroundColor:
            proteinRegionColorMapAnnotations[activeProtein] ||
            "rgba(0, 0, 0, 0.2)",
          borderColor:
            proteinRegionColorMap[activeProtein] || "rgba(0, 0, 0, 0.2)",
          borderWidth: 2,
          label: {
            content: activeProtein,
            enabled: true,
            position: "start",
          },
          z: 10,
        },
      ];
    } else if (showFullAnnotation) {
      return Object.keys(proteinRegions).map((key) => {
        const range = proteinRegions[key];
        const [start, end] = range
          .split("-")
          .map((e) => Math.floor(e / decimateFactor));
        return {
          display: true,
          type: "box",
          xMin: start,
          xMax: end,
          yMin: 0,
          yMax: 3,
          backgroundColor: proteinRegionColorMapAnnotations[key],
          borderColor: proteinRegionColorMap[key],
          borderWidth: 2,
          label: {
            content: key,
            enabled: true,
            position: "start",
          },
          z: 10,
        };
      });
    }
    return [];
  };

  console.log("current decimatefactor is this: ", currentDecimateFactor);

  useEffect(() => {
    const ctx = chartRef.current?.getContext("2d");
    if (!ctx) return;

    // Destroy the previous chart instance if it exists
    if (chartRef.current.chartInstance) {
      chartRef.current.chartInstance.destroy();
    }

    let startPosition = 0;
    let endPosition = 30000;
    // Define your zoom threshold

    // let currentDecimateFactor = selectedProteinRegion ? 1 : decimateFactor;

    if (selectedProteinRegion && proteinRegions[selectedProteinRegion]) {
      [startPosition, endPosition] = proteinRegions[selectedProteinRegion]
        .split("-")
        .map(Number);
    }

    const decimatedLabels = Array.from(
      {
        length: Math.ceil(
          (endPosition - startPosition) / currentDecimateFactor
        ),
      },
      (_, idx) => {
        const position = startPosition + idx * currentDecimateFactor;
        const nucleotide = genomeSequence[position] || "N";
        return `${position}-${nucleotide}`;
      }
    );

    const noMutationData = [];
    const mutationData = Array(4)
      .fill()
      .map(() => []);

    decimatedLabels.forEach((_, idx) => {
      const position = startPosition + 1 + idx * currentDecimateFactor;
      const refNucleotide = genomeSequence[position] || "N";
      const refIndex = nucleotides.indexOf(refNucleotide);

      const positionProbs = [0, 0, 0, 0];
      const dataIndex = Math.floor(idx / decimateFactor);

      nucleotides.forEach((_, nucIndex) => {
        if (dataIndex < decimatedData[nucIndex].length) {
          positionProbs[nucIndex] = decimatedData[nucIndex][dataIndex];
        }
      });

      noMutationData.push(positionProbs[refIndex]);

      nucleotides.forEach((_, nucIndex) => {
        if (nucIndex === refIndex) {
          mutationData[nucIndex].push(0);
        } else {
          mutationData[nucIndex].push(positionProbs[nucIndex]);
        }
      });
    });

    // console.log("genomeData before fetch:", genomeData);
    // const fetchData = (min, max) => {
    //   if (!genomeData || genomeData.length === 0) {
    //     console.log("Empty genomeData detected");
    //     return []; // Handle empty data case
    //   }

    //   return genomeData.map((dataset, idx) => {
    //     const slicedData = dataset.slice(min, max);
    //     console.log(`Dataset ${idx}:`, slicedData);
    //     if (!slicedData.length) {
    //       return { label: nucleotides[idx], data: Array(max - min).fill(0) };
    //     }
    //     return { label: nucleotides[idx], data: slicedData };
    //   });
    // };

    const data = {
      labels: decimatedLabels,
      datasets: [
        {
          label: "No Mutation",
          data: noMutationData,
          backgroundColor: "rgba(128, 128, 128, 0.8)",
          borderColor: "rgba(128, 128, 128, 1)",
          maxBarThickness: 30,
        },
        ...nucleotides.map((nuc, idx) => ({
          label: `${nuc}`,
          data: mutationData[idx],
          backgroundColor: getColorForNucleotide(nuc),
          borderColor: getColorForNucleotide(nuc),
          maxBarThickness: 30,
        })),
      ],
    };

    const options = {
      animation: false,
      responsive: true,
      maintainAspectRatio: false,
      scales: {
        x: {
          stacked: true,
          min: 0, // Always start from the beginning of the genome
          max: 30000,
          bar: {
            categoryPercentage: 1.0,
            barPercentage: 1.0,
          },
        },
        y: {
          stacked: true,
          max: 1,
          ticks: { stepSize: 0.5 },
        },
      },
      plugins: {
        tooltip: {
          callbacks: {
            title: function (tooltipItems) {
              if (!tooltipItems?.length) return ""; // No tooltip items at all

              const tooltipItem = tooltipItems[0];
              const chart = tooltipItem.chart;

              // Safely get the label for the hovered data point
              const label = chart.data.labels?.[tooltipItem.dataIndex];
              // If label is missing/undefined/null, return empty to avoid errors
              if (!label) return "";

              return label; // or a fallback like "Unknown position" if you prefer
            },
            label: function (tooltipItem) {
              const chart = tooltipItem.chart;

              // Safely get the label
              const positionLabel = chart.data.labels?.[tooltipItem.dataIndex];
              if (!positionLabel) {
                // If we don't have a valid label, return empty to avoid the error
                return "";
              }

              // Extract position from the label (e.g., "1000-A")
              const positionStr = positionLabel.split("-")[0];
              const position = parseInt(positionStr, 10);
              // If parse fails or position is NaN, just return empty
              if (isNaN(position)) {
                return "";
              }

              const refNucleotide = genomeSequence?.[position] || "N";

              const datasetLabel = tooltipItem.dataset?.label;
              const value = tooltipItem.raw;

              // If for some reason raw is missing, just return empty
              if (value == null) return "";

              // Now build the tooltip text
              if (datasetLabel === "No Mutation") {
                return `No mutation (${refNucleotide} → ${refNucleotide}): ${value.toFixed(
                  3
                )}`;
              } else {
                const mutatedNucleotide = datasetLabel || "Unknown";
                return `Mutation (${refNucleotide} → ${mutatedNucleotide}): ${value.toFixed(
                  3
                )}`;
              }
            },
          },
        },

        zoom: {
          zoom: {
            wheel: { enabled: true },
            pinch: { enabled: true },
            mode: "x",
            onZoomComplete: ({ chart }) => {
              try {
                if (!chart || !chart.scales || !chart.scales.x) {
                  console.error("Chart or scales are not ready yet");
                  return;
                }
                const scales = chart.scales?.x || {};
                const min = scales.min ?? 1;
                const max = scales.max ?? 30000;
                const zoomLevel = Math.round(chart.getZoomLevel());
                console.log("scales", scales);
                if (isNaN(zoomLevel) || zoomLevel < 1) {
                  console.error("Invalid zoom level:", zoomLevel);
                  return; // Avoid further processing on invalid zoom levels
                }

                setZoomLevel(zoomLevel);

                console.log(`ZOOM-LEVEL: ${zoomLevel}`);

                // console.log(
                //   `Zoom Range : ${zoomLevel}, Min: ${min} - Max- ${max}`
                // );
                let startPos, endPos;
                if (
                  selectedProteinRegion &&
                  proteinRegions[selectedProteinRegion]
                ) {
                  [startPos, endPos] = proteinRegions[selectedProteinRegion]
                    .split("-")
                    .map(Number);
                } else {
                  startPos = Math.max(
                    0,
                    Math.floor(min * currentDecimateFactor)
                  );
                  endPos = Math.min(
                    genomeSequence.length,
                    Math.floor((max + 100) * currentDecimateFactor)
                  );
                }
                startPos = Math.max(0, Math.floor(startPos));
                endPos = Math.min(genomeSequence.length, Math.ceil(endPos));
                // Function to generate chart datasets to avoid repetition
                const generateChartDatasets = (
                  noMutationData,
                  mutationData,
                  nucleotides
                ) => [
                  {
                    label: "No Mutation",
                    data: noMutationData,
                    backgroundColor: "rgba(128, 128, 128, 0.8)",
                    borderColor: "rgba(128, 128, 128, 1)",
                  },
                  ...nucleotides.map((nuc, idx) => ({
                    label: `${nuc}`,
                    data: mutationData[idx],
                    backgroundColor: getColorForNucleotide(nuc),
                    borderColor: getColorForNucleotide(nuc),
                  })),
                ];

                // Function to calculate mutation data
                const calculateMutationData = (
                  fullResolutionData,
                  nucleotides
                ) =>
                  nucleotides.map((nuc, nucIdx) =>
                    fullResolutionData[nucIdx].map((value) => value ?? 0)
                  );

                // Function to calculate no mutation data
                const calculateNoMutationData = (
                  fullResolutionData,
                  genomeSequence,
                  startPos
                ) =>
                  fullResolutionData[0].map((_, posIdx) => {
                    const refNucleotide =
                      genomeSequence[startPos + posIdx] ?? "N";
                    const refIndex = nucleotides.indexOf(refNucleotide);
                    return fullResolutionData[refIndex][posIdx] ?? 0;
                  });

                // Main zoom handling logic
                if (
                  zoomLevel >= zoomThreshold &&
                  selectedProteinRegion === null
                ) {
                  chunkView = false;

                  // Update full-resolution view
                  const fullResolutionData = genomeData.map((dataset) =>
                    dataset.slice(startPos, endPos)
                  );

                  const fullResolutionLabels = Array.from(
                    { length: fullResolutionData[0].length },
                    (_, idx) => {
                      const position = Math.floor(
                        scales.min * currentDecimateFactor + idx
                      );
                      const nucleotide = genomeSequence[position] || "N";
                      return `${position}-${nucleotide}`;
                    }
                  );

                  const fullResolutionMutationData = calculateMutationData(
                    fullResolutionData,
                    nucleotides
                  );
                  const fullResolutionNoMutationData = calculateNoMutationData(
                    fullResolutionData,
                    genomeSequence,
                    startPos
                  );

                  chart.data.labels = fullResolutionLabels;
                  chart.data.datasets = generateChartDatasets(
                    fullResolutionNoMutationData,
                    fullResolutionMutationData,
                    nucleotides
                  );

                  chart.update();
                } else if (
                  zoomLevel >= zoomThreshold &&
                  proteinRegions[selectedProteinRegion] &&
                  selectedProteinRegion === "ORF1ab"
                ) {
                  chunkView = false;

                  // Handle ORF1ab specific view
                  const fullResolutionData = genomeData.map((dataset) =>
                    dataset.slice(startPos, endPos)
                  );

                  const fullResolutionLabels = Array.from(
                    { length: fullResolutionData[0].length },
                    (_, idx) => {
                      const position = Math.floor(
                        scales.min * currentDecimateFactor + idx
                      );
                      const nucleotide = genomeSequence[position] || "N";
                      return `${position}-${nucleotide}`;
                    }
                  );

                  const fullResolutionMutationData = calculateMutationData(
                    fullResolutionData,
                    nucleotides
                  );
                  const fullResolutionNoMutationData = calculateNoMutationData(
                    fullResolutionData,
                    genomeSequence,
                    startPos
                  );

                  chart.data.labels = fullResolutionLabels;
                  chart.data.datasets = generateChartDatasets(
                    fullResolutionNoMutationData,
                    fullResolutionMutationData,
                    nucleotides
                  );

                  chart.update();
                }

                // Single update call outside the conditions
                else {
                  chunkView = true;

                  // Update decimated view
                  chart.data.labels = decimatedLabels;
                  chart.data.datasets = [
                    {
                      label: "No Mutation",
                      data: noMutationData,
                      backgroundColor: "rgba(128, 128, 128, 0.8)",
                      borderColor: "rgba(128, 128, 128, 1)",
                    },
                    ...nucleotides.map((nuc, idx) => ({
                      label: `${nuc}`,
                      data: mutationData[idx],
                      backgroundColor: getColorForNucleotide(nuc),
                      borderColor: getColorForNucleotide(nuc),
                    })),
                  ];

                  chart.update();
                }
              } catch (error) {
                console.error("Error during zoom handling:", error);
                console.error("Detailed Zoom Handling Error:", error);
                console.error("Error Name:", error.name);
                console.error("Error Message:", error.message);
                console.error("Error Stack:", error.stack);
              }
            },
          },
          pan: {
            enabled: true,
            mode: "x",
            threshold: 15,
          },
          limits: {
            x: { min: 30, max: 300000 + 100, minRange: 25, maxRange: 65000 },
          },
        },

        annotation: {
          annotations: createAnnotations(),
        },
      },
    };

    const chartInstance = new Chart(ctx, {
      type: "bar",
      data,
      options,
    });

    chartRef.current.chartInstance = chartInstance;

    chartInstance.resetZoom();

    return () => {
      chartInstance.destroy();
    };
  }, [decimatedData, genomeSequence, selectedProteinRegion]);

  useEffect(() => {
    const chartInstance = chartRef.current.chartInstance;
    if (chartInstance) {
      const annotations = createAnnotations();
      chartInstance.options.plugins.annotation.annotations = annotations;
      chartInstance.update();
    }
  }, [activeProtein, showFullAnnotation, maxValue]);

  const handleProteinHover = (protein) => setActiveProtein(protein);
  const handleProteinLeave = () => setActiveProtein(null);
  const handleShowFullAnnotation = () => setShowFullAnnotation((prev) => !prev);

  return (
    <div className="overflow-x-hidden">
      <div className="chart-container w-full flex bg-[#f6f7f9] py-5">
        <SidePanel
          proteinRegions={proteinRegions}
          onProteinHover={handleProteinHover}
          onProteinLeave={handleProteinLeave}
          handleShowFullAnnotation={handleShowFullAnnotation}
        />
        <div className="w-full">
          <canvas
            className="w-full h-[90vh] max-h-screen bg-white mt-6 pl-4 pr-8 py-2 rounded-xl shadow-md"
            ref={chartRef}
          />
        </div>
      </div>
    </div>
  );
};

export default GenomeChart;
