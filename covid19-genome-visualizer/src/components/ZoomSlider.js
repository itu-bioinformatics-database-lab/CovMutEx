import React, { useState } from "react";
import { FiPlusCircle, FiMinusCircle } from "react-icons/fi";

const ZoomSlider = ({ handleZoom }) => {
	const [zoomLevel, setZoomLevel] = useState(10); // Adjust initial zoom level as needed

	const handleZoomChange = (event) => {
		setZoomLevel(event.target.value);

		handleZoom(event.target.value / 10);
	};

	const zoomIn = () => {
		setZoomLevel((zoomLevel) => Math.min(zoomLevel + 1, 20));
		handleZoom(1.1);
	};

	const zoomOut = () => {
		setZoomLevel((zoomLevel) => Math.max(zoomLevel - 1, 0));
		handleZoom(0.9);
	};

	return (
		<div className="zoom-container">
			<button onClick={zoomIn}>
				<FiPlusCircle color="#333333" className="w-8 h-8" />
			</button>
			<input
				type="range"
				className="zoom-slider"
				min="0"
				max="20"
				value={zoomLevel}
				onChange={handleZoomChange}
				orient="vertical"
			/>
			<button onClick={zoomOut}>
				<FiMinusCircle color="#333333" className="w-8 h-8" />
			</button>
		</div>
	);
};

export default ZoomSlider;
