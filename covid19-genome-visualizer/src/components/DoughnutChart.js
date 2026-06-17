import React, { useMemo, useState } from "react";
import { Doughnut } from "react-chartjs-2";
import { useDispatch, useSelector } from "react-redux";
import { showProteinRegion } from "../features/genome/genomeSlice";
import { proteinRegionColorMap } from "../utils/proteinRegionColorMap";
import { proteinRegionsSize as defaultProteinRegionsSize } from "../data/proteinRegions";
import { Switch } from "@material-tailwind/react";

const DoughnutChart = ({ data = {} }) => {
  // All hooks must be called unconditionally and in the same order on every
  // render. Any early return (e.g. the "no data" guard below) MUST come AFTER
  // every hook call here, per the React rules-of-hooks.
  const [normalized, setNormalized] = useState(false);
  const [title, setTitle] = useState("Mutation Probability");
  const dispatch = useDispatch();

  // Organism-aware protein region lengths. proteinRegionPossibilities is
  // {name: [start, end]} (1-based inclusive) coming from the backend's
  // organism layer; length = end - start + 1. Falls back to the hardcoded
  // COVID sizes when no prediction has landed yet.
  const organismRegionPossibilities = useSelector(
    (state) => state.genome.proteinRegionPossibilities
  );
  const proteinRegionsSize = useMemo(() => {
    if (
      organismRegionPossibilities &&
      Object.keys(organismRegionPossibilities).length > 0
    ) {
      return Object.fromEntries(
        Object.entries(organismRegionPossibilities).map(([name, range]) => [
          name,
          Array.isArray(range) ? range[1] - range[0] + 1 : 1,
        ])
      );
    }
    return defaultProteinRegionsSize;
  }, [organismRegionPossibilities]);

  // Don't render if no data — early return must come AFTER all hooks above.
  if (!data || Object.keys(data).length === 0) {
    return null;
  }

  function calculateTotalSum() {
    return Object.entries(data).reduce((sum, [key, value]) => {
      const length = proteinRegionsSize[key] || 1;
      const normalizedValue = normalized ? value / length : value;
      return sum + normalizedValue;
    }, 0);
  }
  const totalSum = calculateTotalSum();
  const percentages = {};
  Object.keys(data).forEach((key) => {
    const length = proteinRegionsSize[key] || 1;
    const value = data[key] || 0;
    const normalizedValue = normalized ? value / length : value;
    const percentage = totalSum > 0 ? (normalizedValue / totalSum) * 100 : 0;
    percentages[key] = parseFloat(percentage.toFixed(2));
  });

  const handleNormalizedButton = () => {
    setTitle(
      normalized ? "Mutation Probability" : "Mutation Probability (Normalized)"
    );
    setNormalized(!normalized);
  };

  const chartData = {
    labels: Object.keys(percentages),
    datasets: [
      {
        label: "Mutation Probability (%)",
        data: Object.values(percentages),
        backgroundColor: Object.keys(percentages).map(
          (key) => proteinRegionColorMap[key] || "#ccc"
        ),
        borderWidth: 1,
      },
    ],
  };

  const chartOptions = {
    responsive: true,
    maintainAspectRatio: true,
    plugins: {
      legend: { position: "top", labels: { boxWidth: 12, padding: 8, font: { size: 11 } } },
      title: { display: false },
    },
    onClick: (event, elements) => {
      if (elements[0]) {
        const clickedIndex = elements[0].index;
        const clickedLabel = chartData.labels[clickedIndex];
        // dispatch(showProteinRegion(clickedLabel));
      }
    },
  };

  return (
    <div className="w-full max-w-5xl mx-auto bg-white dark:bg-gray-900 border border-gray-200 dark:border-gray-700 shadow-md rounded-xl p-6">
      {/* Header with normalize toggle */}
      <div className="flex items-center justify-between mb-4">
        <h2 className="text-lg font-bold text-gray-800 dark:text-gray-200">
          Protein Mutation Probability
        </h2>
        <div className="flex items-center gap-2 cursor-pointer">
          <span className="text-sm font-medium text-gray-600 dark:text-gray-400">
            Normalize by region length
          </span>
          <Switch
            onClick={handleNormalizedButton}
            containerProps={{ className: "!bg-gray-300" }}
            color="blue"
            className="custom-switch"
          />
        </div>
      </div>

      {/* Horizontal layout: Doughnut + Table */}
      <div className="flex flex-col md:flex-row items-start gap-6">
        {/* Doughnut Chart */}
        <div className="w-full md:w-1/2 flex justify-center">
          <div className="w-[320px] max-w-full">
            <Doughnut data={chartData} options={chartOptions} />
          </div>
        </div>

        {/* Data Table */}
        <div className="w-full md:w-1/2">
          <h3 className="text-sm font-semibold text-gray-800 dark:text-gray-200 mb-2">
            Data View
          </h3>
          <div className="overflow-hidden rounded-lg border border-gray-300 dark:border-gray-600">
            <table className="w-full text-sm text-left text-gray-600 dark:text-gray-400">
              <thead className="text-xs text-gray-700 dark:text-gray-300 uppercase bg-gray-50 dark:bg-gray-800">
                <tr>
                  <th scope="col" className="px-3 py-2 border-r border-gray-300 dark:border-gray-600">
                    Protein Region
                  </th>
                  <th scope="col" className="px-3 py-2 text-right">
                    Prob. (%)
                  </th>
                </tr>
              </thead>
              <tbody>
                {Object.entries(percentages).map(([region, percentage]) => (
                  <tr
                    key={region}
                    className="bg-white dark:bg-gray-900 border-b border-gray-200 dark:border-gray-700 hover:bg-gray-50 dark:hover:bg-gray-800 transition-colors"
                  >
                    <td className="px-3 py-1.5 font-medium text-gray-800 dark:text-gray-200 whitespace-nowrap border-r border-gray-200 dark:border-gray-700">
                      <div className="flex items-center gap-2">
                        <span
                          className="w-3 h-3 rounded-full inline-block flex-shrink-0"
                          style={{ backgroundColor: proteinRegionColorMap[region] || "#ccc" }}
                        />
                        {region}
                      </div>
                    </td>
                    <td className="px-3 py-1.5 text-right font-mono">
                      {percentage}%
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      </div>
    </div>
  );
};

export default DoughnutChart;
