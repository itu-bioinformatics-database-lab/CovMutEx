import { Progress } from "@material-tailwind/react";
import React, { useState, useEffect } from "react";

export default function Progressbar() {
	const [filled, setFilled] = useState(0);
	const [isRunning, setIsRunning] = useState(true);
	useEffect(() => {
		if (filled < 100 && isRunning) {
			if (filled < 96) {
				setTimeout(() => setFilled((prev) => (prev += 1)), 400);
			}
		}
	}, [filled, isRunning]);

	useEffect(() => {
		if (filled >= 100) {
			setIsRunning(false);
		}
	}, [filled]);
	return (
		<div className="grid h-screen place-items-center">
			<div className="w-[400px] flex flex-col items-center">
				<p className="mt-6 text-md leading-8 text-gray-900 italic">
					Predicting...
				</p>
				<Progress size="sm" value={filled} color="blue" />
			</div>
		</div>
	);
}
