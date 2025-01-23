import React, { useRef, useEffect } from "react";
import Chart from "chart.js/auto";
import "chartjs-plugin-zoom";

const ZoomChart = () => {
	const chartContainer = useRef(null);
	let timer;

	useEffect(() => {
		const zoomOptions = {
			pan: {
				enabled: true,
				mode: "x",
				modifierKey: "ctrl",
				onPanComplete: startFetch,
			},
			zoom: {
				wheel: {
					enabled: true,
				},
				drag: {
					enabled: true,
				},
				pinch: {
					enabled: true,
				},
				mode: "x",
				onZoomComplete: startFetch,
			},
		};

		const allData = Array.from({ length: 30000 }, (_, index) => ({
			x: index,
			y: Math.random(),
		}));

		function fetchData(x1, x2) {
			return allData.slice(x1, x2);
		}

		function startFetch({ chart }) {
			let { min, max } = chart.scales.x;
			min = parseInt(min);
			max = parseInt(max);
			clearTimeout(timer);
			timer = setTimeout(() => {
				// console.log("Fetched data between " + min + " and " + max);
				chart.data.datasets[0].data = fetchData(min, max);
				chart.stop(); // make sure animations are not running
				chart.update("none");
			}, 500);
		}

		const config = {
			type: "bar",
			data: {
				datasets: [
					{
						label: "Random Data",
						borderColor: "red",
						backgroundColor: "rgba(255, 0, 0, 0.1)",
						data: fetchData(0, 500),
					},
				],
			},
			options: {
				scales: {
					x: {
						position: "bottom",
						type: "linear",
						ticks: {
							autoSkip: true,
							autoSkipPadding: 50,
							maxRotation: 0,
							callback: function (value) {
								// Format the label
								return `Position ${parseInt(value)}`;
							},
						},
					},
					y: {
						type: "linear",
						position: "left",
						min: 0,
						max: 1,
					},
				},
				plugins: {
					zoom: zoomOptions,
					title: {
						display: true,
						position: "bottom",
						text: "Zoomable Bar Chart",
					},
				},
			},
		};

		const ctx = chartContainer.current.getContext("2d");
		new Chart(ctx, config);
	}, []);

	return (
		<div className="w-[100vw] h-[80vh] mx-5">
			<canvas ref={chartContainer} width="800" height="400"></canvas>
		</div>
	);
};

export default ZoomChart;
