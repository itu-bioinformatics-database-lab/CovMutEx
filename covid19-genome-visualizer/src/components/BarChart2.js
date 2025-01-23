import React, { useEffect, useRef, useState } from "react";
import Chart from "chart.js/auto";
import zoomPlugin from "chartjs-plugin-zoom";
import { useDispatch, useSelector } from "react-redux";
import { resetChart } from "../features/genome/genomeSlice";
import { proteinRegionColorMap } from "../utils/proteinRegionColorMap";
import { GrPowerReset } from "react-icons/gr";
import { Button } from "@material-tailwind/react";
import ZoomSlider from "./ZoomSlider";
Chart.register(zoomPlugin);
let nucleotides = { 0: "A", 1: "C", 2: "T", 3: "G" };

function BarChart2({ data, seq }) {
	const chartRef = useRef();
	const dispatch = useDispatch();
	const chartTitle = useSelector((state) => state.genome.chartTitle);
	const isWholeSequenceSelected = useSelector(
		(state) => state.genome.isWholeSequenceSelected
	);

	const [_chart, setChart] = useState(null);
	const [click, setClick] = useState(false);
	// const decimateData = (data, factor) => {
	// 	console.log("dtaaaa", data);
	// 	return data.reduce((acc, _, index) => {
	// 		if (index % factor === 0) {
	// 			const slice = data.slice(index, index + factor);
	// 			const averaged = slice.reduce((avg, curr) => {
	// 				Object.keys(curr.mutationPoss).forEach((key) => {
	// 					avg[key] = (avg[key] || 0) + curr.mutationPoss[key] / factor;
	// 				});
	// 				return avg;
	// 			}, {});
	// 			acc.push({
	// 				pos: index / factor,
	// 				nucleotide: "",
	// 				mutationPoss: averaged,
	// 			});
	// 		}
	// 		return acc;
	// 	}, []);
	// };
	const decimateData = (data, factor) => {
		return data.map((dataset) =>
			dataset.reduce((acc, curr, index) => {
				if (index % factor === 0) {
					const slice = dataset.slice(index, index + factor);
					const sum = slice.reduce((sum, value) => sum + value, 0);
					acc.push(sum / factor);
				}
				return acc;
			}, [])
		);
	};

	useEffect(() => {
		const decimatedData = decimateData(data, 15); // Adjust factor as needed
		const labels = seq.slice(0, decimatedData[0].length);
		const label_indexes = labels.split("").map((_, idx) => idx);
		console.log("lvlll", label_indexes);
		// for (const row of data) {
		// 	dataset.push(row.slice(0, 1000).map((num) => parseFloat(num)));
		// }
		// console.log("dataset:", dataset);
		const currentChartRef = chartRef.current;
		// const decimatedData = decimateData(data, 30); // Adjust factor as needed

		let idx = -1;
		const datasets = decimatedData.map((dataset, idx) => ({
			label: nucleotides[idx],
			data: dataset,
			borderColor: getColorForNucleotide(nucleotides[idx]),
			backgroundColor: getColorForNucleotide(nucleotides[idx]),
		}));
		console.log("datasets:", datasets);
		console.log("labels", labels);

		if (currentChartRef) {
			const ctx = currentChartRef.getContext("2d");

			if (ctx) {
				const chart = new Chart(ctx, {
					type: "bar",
					data: {
						labels: label_indexes,
						datasets: datasets,
					},
					options: {
						animation: false,
						scales: {
							x: {
								// id: "xAxis1",
								// type: "category",
								stacked: true,
								// labels: labels,
							},

							y: {
								// beginAtZero: true,
								stacked: true,
								// max: 4,
							},
						},
						// tooltips: {
						// 	mode: "index",
						// 	intersect: false,
						// },
						// hover: {
						// 	mode: "index",
						// 	intersect: false,
						// },
						responsive: true,
						maintainAspectRatio: false,
						plugins: {
							zoom: {
								zoom: {
									wheel: {
										enabled: true,
									},
									pinch: {
										enabled: true,
									},
									mode: "x",
								},
								pan: {
									enabled: true,
									mode: "x",
									onPan: () => {
										console.log("panned");
									},
								},
							},
							title: {
								display: true,
								color: proteinRegionColorMap[chartTitle],
								position: "bottom",
								text: chartTitle,
								font: {
									size: 16,
								},
							},
						},
					},
				});
				setChart(chart);

				return () => chart.destroy();
			}
		}
	}, [data, seq]);

	const handleReset = () => {
		dispatch(resetChart());
	};

	const handleZoom = (zoomLevel) => {
		console.log("zoom:", zoomLevel);
		_chart.zoom(zoomLevel);
	};

	return (
		<div className="flex">
			<div className="flex items-center justify-center chart-container">
				<ZoomSlider handleZoom={handleZoom} />
			</div>

			<div
				className="chart-container"
				style={{ height: "50vh", width: "60vw" }}
			>
				<div className=" ml-12">
					<Button
						color="blue"
						variant="gradient"
						type="submit"
						className="flex gap-2 justify-center items-center"
						onClick={handleReset}
					>
						<div>Reset Chart</div>
						<GrPowerReset size={20} color="white" />
					</Button>
				</div>

				<canvas ref={chartRef} />
				{/* <div className="flex">
				<button
					className="border-4 border-gray-400 m-2 p-2"
					onClick={(e) => {
						console.log(_chart);
						_chart.zoom(1.1);
					}}
				>
					{" "}
					+
				</button>
				<button
					className="border-4 border-gray-400 m-2 p-2"
					onClick={(e) => {
						console.log(_chart);
						_chart.zoom(0.9);
					}}
				>
					{" "}
					-
				</button>
			</div> */}
			</div>
		</div>
	);
}

function getColorForNucleotide(nucleotide) {
	const colorMap = {
		A: "#FF5733",
		C: "#3399FF",
		T: "#33CC33",
		G: "#9966FF",
	};

	return colorMap[nucleotide] || "#000000";
}

export default BarChart2;
