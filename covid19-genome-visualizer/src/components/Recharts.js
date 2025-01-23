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

const GenomeChart = ({ genomeData, genomeSequence }) => {
  const chartRef = useRef(null);
  const [activeProtein, setActiveProtein] = useState(null);
  const [showFullAnnotation, setShowFullAnnotation] = useState(false);
  const [hoveredPosition, setHoveredPosition] = useState(null);
  const selectedProteinRegion = useSelector(
    (state) => state.genome.selectedProteinRegion
  );
  const genomeState = useSelector((state) => state.genome);
  const decimateFactor = 30;

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

  useEffect(() => {
    const ctx = chartRef.current?.getContext("2d");
    if (!ctx) return;

    // Destroy the previous chart instance if it exists
    if (chartRef.current.chartInstance) {
      chartRef.current.chartInstance.destroy();
    }

    let startPosition = 0;
    let endPosition = 30000;

    if (selectedProteinRegion && proteinRegions[selectedProteinRegion]) {
      [startPosition, endPosition] = proteinRegions[selectedProteinRegion]
        .split("-")
        .map(Number);
    }

    const decimatedLabels = Array.from(
      { length: Math.ceil((endPosition - startPosition) / decimateFactor) },
      (_, idx) => {
        const position = startPosition + idx * decimateFactor;
        const nucleotide = genomeSequence[position] || "N";
        return `${position}-${nucleotide}`;
      }
    );

    const noMutationData = [];
    const mutationData = Array(4)
      .fill()
      .map(() => []);

    decimatedLabels.forEach((_, idx) => {
      const position = startPosition + idx * decimateFactor;
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

    const data = {
      labels: decimatedLabels,
      datasets: [
        {
          label: "No Mutation",
          data: noMutationData,
          backgroundColor: "rgba(128, 128, 128, 0.8)",
          borderColor: "rgba(128, 128, 128, 1)",
          maxBarThickness: 20,
        },
        ...nucleotides.map((nuc, idx) => ({
          label: `${nuc}`,
          data: mutationData[idx],
          backgroundColor: getColorForNucleotide(nuc),
          borderColor: getColorForNucleotide(nuc),
          maxBarThickness: 20,
        })),
      ],
    };

    const options = {
      animation: true,
      responsive: true,
      maintainAspectRatio: false,
      scales: {
        x: {
          stacked: true,
          min: 0, // Always start from the beginning of the genome
          max: 30000, // Always end at the full genome range
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
              const idx = tooltipItems[0].dataIndex;
              const position = startPosition + idx * decimateFactor;
              const nucleotide = genomeSequence[position] || "N";

              if (selectedProteinRegion) {
                return `${selectedProteinRegion} Position ${position} (Reference: ${nucleotide})`;
              }
              return `Position ${position} (Reference: ${nucleotide})`;
            },
            label: function (tooltipItem) {
              const idx = tooltipItem.dataIndex;
              const position = startPosition + idx * decimateFactor;
              const refNucleotide = genomeSequence[position] || "N";

              if (tooltipItem.dataset.label === "No Mutation") {
                return `No mutation (${refNucleotide} → ${refNucleotide}): ${tooltipItem.raw.toFixed(
                  3
                )}`;
              } else {
                const mutatedNucleotide = tooltipItem.dataset.label;
                return `Mutation (${refNucleotide} → ${mutatedNucleotide}): ${tooltipItem.raw.toFixed(
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
          },
          pan: {
            enabled: true,
            mode: "x",
          },
          limits: {
            x: { min: 0, max: 30000 },
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

    // Reset zoom explicitly to the default range
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
