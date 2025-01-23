import React, { useEffect, useRef, useState } from "react";
import Chart from "chart.js/auto";
import zoomPlugin from "chartjs-plugin-zoom";
import { useDispatch, useSelector } from "react-redux";
import { resetChart } from "../features/genome/genomeSlice";
import { proteinRegionColorMap } from "../utils/proteinRegionColorMap";
import { proteinRegions } from "../data/proteinRegions";
Chart.register(zoomPlugin);

// function chunkAndAverageData(data, chunkSize) {
// 	const chunkedData = [];
// 	for (let i = 0; i < data.length; i += chunkSize) {
// 		let chunk = data.slice(i, i + chunkSize);
// 		let sum = chunk
// 			.map((item) => item["mutationPoss"])
// 			.reduce((acc, { A, T, G, C }) => acc + A + T + G + C, 0);
// 		let avg = parseFloat((sum / chunkSize).toFixed(2));
// 		chunkedData.push(avg);
// 	}
// 	return chunkedData;
// }

function BarChart({ data }) {
	// console.log("barchartdata: ", data);
	const chartRef = useRef();
	const dispatch = useDispatch();
	const [chunkView, setChunkView] = useState(true);
	const chartTitle = useSelector((state) => state.genome.chartTitle);
	const [_chart, setChart] = useState(null);
	const isWholeSequenceSelected = useSelector(
		(state) => state.genome.isWholeSequenceSelected
	);

	const decimateData = (data, factor) => {
		console.log("dtaaaa", data);
		return data.reduce((acc, _, index) => {
			if (index % factor === 0) {
				const slice = data.slice(index, index + factor);
				const averaged = slice.reduce((avg, curr) => {
					Object.keys(curr.mutationPoss).forEach((key) => {
						avg[key] = (avg[key] || 0) + curr.mutationPoss[key] / factor;
					});
					return avg;
				}, {});
				acc.push({
					pos: index / factor,
					nucleotide: "",
					mutationPoss: averaged,
				});
			}
			return acc;
		}, []);
	};
	useEffect(() => {
		const currentChartRef = chartRef.current;
		const decimatedData = decimateData(data, 30); // Adjust factor as needed

		// const chunkSize = 5;
		// const chunkedData = chunkAndAverageData(data, chunkSize);

		// const chunkedDatasets = [
		// 	{
		// 		label: "Average Data",
		// 		data: chunkedData,
		// 		backgroundColor: "rgba(0, 123, 255, 0.5)", // or any color you prefer
		// 	},
		// ];
		// console.log("chunked data: ", chunkedData);

		const datasets = Object.keys(decimatedData[0].mutationPoss).map(
			(nucleotide) => ({
				label: nucleotide,
				data: decimatedData.map((entry) => entry.mutationPoss[nucleotide]),
				borderColor: getColorForNucleotide(nucleotide),
				backgroundColor: getColorForNucleotide(nucleotide),
				// stack: "stack",
			})
		);

		if (currentChartRef) {
			const ctx = currentChartRef.getContext("2d");

			if (ctx) {
				const labels = decimatedData.map(
					(entry) => `${entry.pos}-${entry.nucleotide}`
				);
				// const labels2 = data.map(
				// 	(entry) => `${entry.pos}-${entry.nucleotide}-${entry.proteinRegion}`
				// );

				// const labels_chunk = chunkedData.map((chunk, idx) => `Chunk-${idx}`);

				// const datasets =

				// console.log("1111:");
				// console.log(labels);
				// console.log(datasets);
				const chart = new Chart(ctx, {
					type: "bar",
					data: {
						labels: labels,
						datasets: datasets,
					},
					options: {
						animation: false,
						scales: {
							x: {
								id: "xAxis1",
								type: "category",
								stacked: true,
								labels: labels,
								// ticks: {
								// 	callback: function (val, index) {
								// 		// Hide the label of every 2nd dataset
								// 		const label = this.getLabelForValue(val);
								// 		return label.split("-")[1];
								// 	},
								// },
							},
							// x2: {
							// 	id: "xAxis2",
							// 	type: "category",
							// 	gridLines: {
							// 		display: false,
							// 	},
							// 	labels: labels2,
							// 	ticks: {
							// 		autoSkip: true,
							// 		// callback: function (val, index) {
							// 		// 	// Hide the label of every 2nd dataset
							// 		// 	const pR = this.getLabelForValue(val).split("-")[2];
							// 		// 	// console.log("pR", pR);
							// 		// 	const index_ = parseInt(
							// 		// 		this.getLabelForValue(val).split("-")[0]
							// 		// 	);
							// 		// 	// console.log(proteinRegions[pR.split(" ")[0].split("-")]);
							// 		// 	const [start, end] =
							// 		// 		proteinRegions[pR.split(" ")[0].split("-")].split("-");
							// 		// 	const middle = parseInt(
							// 		// 		(parseInt(start) + parseInt(end)) / 2
							// 		// 	);
							// 		// 	if (isWholeSequenceSelected) {
							// 		// 		if (index_ % 100 === 0 && index_ !== 0) {
							// 		// 			// console.log("YES");
							// 		// 			return pR;
							// 		// 		}
							// 		// 	} else {
							// 		// 		if (index_ === middle) {
							// 		// 			return pR;
							// 		// 		}
							// 		// 	}
							// 		// },
							// 	},
							// },

							y: {
								beginAtZero: true,
								stacked: true,
								max: 4,
							},
						},
						tooltips: {
							mode: "index",
							intersect: false,
						},
						hover: {
							mode: "index",
							intersect: false,
						},
						responsive: true,
						maintainAspectRatio: true,
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
									// onZoomComplete: (context) => {
									// 	console.log("AAAAAAAAAAAAAAAAA");
									// 	console.log(context.chart.getZoomLevel());
									// 	if (context.chart.getZoomLevel() >= 1.4) {
									// 		setChunkView((_) => false); // Switch to detailed view

									// 		context.chart.data.datasets = detailedDatasets;
									// 		console.log("worked onzoomm");
									// 		context.chart.update();
									// 	} else {
									// 		setChunkView((_) => true);
									// 		context.chart.data.datasets = chunkedDatasets;
									// 		context.chart.update();
									// 	}
									// },
								},
								pan: {
									enabled: true,
									mode: "x",
									onPan: () => {
										console.log("panned");
										// const leftEnd =
										// 	chart.getDatasetMeta(0).dataset._scale.chart.scales["x"]
										// 		._table[0];
										// console.log(chart.getDatasetMeta(0));
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
	}, [data]);

	const handleReset = () => {
		setChunkView(true);
		dispatch(resetChart());
	};

	return (
		<div>
			<button onClick={handleReset}>Reset</button>
			<canvas ref={chartRef} />
			<div className="flex">
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

export default BarChart;
