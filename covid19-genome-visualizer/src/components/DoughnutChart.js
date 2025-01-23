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

  // Calculate the total sum of mutation probabilities
  function calculateTotalSum() {
    return Object.entries(data).reduce((sum, [key, value]) => {
      const length = proteinRegionsSize[key] || 1; // Default to 1 to avoid division by zero
      const normalizedValue = normalized ? value / length : value;
      return sum + normalizedValue;
    }, 0);
  }

  const totalSum = calculateTotalSum(); // Total sum for percentage calculation

  // Calculate percentages for each protein region
  const percentages = {};
  Object.keys(data).forEach((key) => {
    const length = proteinRegionsSize[key] || 1; // Avoid division by zero
    const value = data[key] || 0; // Default to 0 if undefined

    const normalizedValue = normalized ? value / length : value;
    percentages[key] = parseFloat(
      ((normalizedValue / totalSum) * 100).toFixed(2)
    );
  });

  // Handle normalization toggle
  const handleNormalizedButton = () => {
    setTitle(
      normalized ? "Mutation Probability" : "Mutation Probability (Normalized)"
    );
    setNormalized(!normalized);
  };

  // Prepare data for the Doughnut chart
  const chartData = {
    labels: Object.keys(percentages),
    datasets: [
      {
        label: "Mutation Probability (%)",
        data: Object.values(percentages),
        backgroundColor: Object.keys(percentages).map(
          (key) => proteinRegionColorMap[key] || "#ccc" // Default to gray if color is missing
        ),
        borderWidth: 1,
      },
    ],
  };

  // Doughnut chart options
  const chartOptions = {
    responsive: true,
    maintainAspectRatio: true,
    plugins: {
      legend: {
        position: "top",
      },
      title: {
        display: true,
        text: title,
      },
    },
    onClick: (event, elements) => {
      if (elements[0]) {
        const clickedIndex = elements[0].index;
        const clickedLabel = chartData.labels[clickedIndex];
        // dispatch(showProteinRegion(clickedLabel)); // Dispatch to Redux
      }
    },
  };

  return (
    <div className="flex justify-center items-center h-full">
      <div
        className="w-[400px] mt-12 h-[90vh]
                   flex flex-col justify-center
                   text-center bg-white 
                   border-[2px] gap-x-0
                   shadow-md rounded-md py-10"
      >
        <div className="flex justify-start items-center cursor-pointer">
          <div className="ms-3 font-bold flex text-gray-700">Normalize</div>
          <div className="pl-2 pt-1">
            <Switch
              onClick={handleNormalizedButton}
              color="blue"
              className="custom-switch"
            />
          </div>
        </div>
        <Doughnut className="" data={chartData} options={chartOptions} />
      </div>
    </div>
  );
};

export default DoughnutChart;
