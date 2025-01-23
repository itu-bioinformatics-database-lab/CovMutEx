import { Spinner } from "@material-tailwind/react";
import React from "react";

const LoadingSpinner = () => {
	return (
		<div className="grid h-screen place-items-center">
			<Spinner color="blue" className="h-10 w-10" />
		</div>
	);
};

export default LoadingSpinner;
