import React from "react";
import { useSelector } from "react-redux";

const RandomSequence = () => {
	const randomSeq = useSelector((state) => state.genome.randomSeq);

	return (
		<div>
			<h2>Random Genome:</h2>
			<pre>{JSON.stringify(randomSeq)}</pre>
		</div>
	);
};

export default RandomSequence;
