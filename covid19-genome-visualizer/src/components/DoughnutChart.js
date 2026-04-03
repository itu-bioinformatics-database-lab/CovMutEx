import React, { useState } from "react";
import { Doughnut } from "react-chartjs-2";
import { useDispatch } from "react-redux";
import { showProteinRegion } from "../features/genome/genomeSlice";
import { proteinRegionColorMap } from "../utils/proteinRegionColorMap";
import { proteinRegionsSize } from "../data/proteinRegions";
import { Switch } from "@material-tailwind/react";

const DoughnutChart = ({ data }) => {
  const [normalized, setNormalized] = useState(false);
  const [title, setTitle] = useState("Mutation Probability");
  const dispatch = useDispatch();

  // ... (all the calculation logic remains the same)
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
      legend: { position: "top" },
      title: { display: true, text: title },
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
    <div className="flex justify-center items-start h-full">
      {/* --- Main container with fixed height and scroll --- */}
      <div
        className="w-[350px] ms-2 mt-12 h-[90vh] 
                   flex flex-col
                   text-center bg-white 
                   border-[2px]
                   shadow-md rounded-md py-6     
                   overflow-hidden            
                   "
      >
        {/* --- SWITCH AND CHART SECTION --- */}
        <div className="flex-shrink-0">
          <div className="flex justify-center items-center cursor-pointer mb-2">
            <div className="ms-3 font-bold flex text-gray-700">
              Normalize (by protein region length)
            </div>
            <div className="pl-2 pt-1">
              <Switch
                onClick={handleNormalizedButton}
                containerProps={{ className: "!bg-gray-300" }}
                color="blue"
                className="custom-switch"
              />
            </div>
          </div>
          {/* --- DOUGHNUT CHART --- */}
          <div className="px-4">
            <Doughnut data={chartData} options={chartOptions} />
          </div>
        </div>

        {/* --- DATA TABLE SECTION --- */}
        <div className="mt-4 px-4">
          <div className="border-t pt-2">
            <h3 className="text-sm font-semibold text-gray-800 mb-1">
              Data View
            </h3>
            <div>
              <table className="w-full text-xs text-left text-gray-600 border border-gray-400">
                <thead className="text-xs text-gray-700 uppercase bg-gray-50">
                  <tr>
                    <th scope="col" className="px-1 py-1 border-r border-gray-400">
                      Protein Region
                    </th>
                    <th scope="col" className="px-1 py-1 text-right">
                      Prob. (%)
                    </th>
                  </tr>
                </thead>
                <tbody>
                  {Object.entries(percentages)
                    .map(([region, percentage]) => (
                    <tr
                      key={region}
                      className="bg-white border-b border-gray-400 hover:bg-gray-50"
                    >
                      <td className="px-1 py-0.5 font-medium text-gray-800 whitespace-nowrap border-r border-gray-400">
                        <span>{region}</span>
                      </td>
                      <td className="px-1 py-0.5 text-right">
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
    </div>
  );
};

export default DoughnutChart;
