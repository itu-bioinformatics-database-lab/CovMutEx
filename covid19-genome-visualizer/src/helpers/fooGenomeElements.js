import {
	nucleotides,
	windowSize,
	sequenceLength,
} from "../constants/constantVars";
import { proteinRegions } from "../data/proteinRegions";

// parameters: None
// objective: generates 30k long random strain
// return: 30k length long array
export function generateRandomSequence(length = sequenceLength) {
	const randomSequence = Array.from(
		{ length },
		() => nucleotides[Math.floor(Math.random() * nucleotides.length)]
	);
	return randomSequence;
}

export function featureExtractor(kmers = null) {
	return kmers;
}

export function fooModelOutputGenerator(kmerFeatures = null) {
	const randomNumber = Math.random();
	return parseFloat(randomNumber.toFixed(2));
}

// this function is hardcoded for windowSize=5 for now...
// TODO: refactor this function
function windowPropagation(randomGenome, windowLength) {
	const windowSliceArr = [];
	const possibilityMap = randomGenome.map((nucleotide, index) => {
		let windowSlice = [];
		if (index === 0) {
			windowSlice = randomGenome.slice(index, index + windowLength - 2);
			windowSlice.unshift("-");
			windowSlice.unshift("-");
		} else if (index === 1) {
			windowSlice = randomGenome.slice(index - 1, index + windowLength - 2);
			windowSlice.unshift("-");
		} else if (index === 2) {
			windowSlice = randomGenome.slice(index - 2, index + windowLength - 2);
		} else if (index === randomGenome.length - 2) {
			windowSlice = randomGenome.slice(index - 2, index + windowLength - 2);
			windowSlice.push("-");
		} else if (index === randomGenome.length - 1) {
			windowSlice = randomGenome.slice(index - 2, index + windowLength - 2);
			windowSlice.push("-");
			windowSlice.push("-");
		} else {
			windowSlice = randomGenome.slice(index - 2, index + windowLength - 2);
		}
		windowSliceArr.push(windowSlice);
		// windowSlice feature extract
		const mutationPoss = calculateMutationPoss(windowSlice);

		let proteinRegion = determineProteinRegion(index);

		return {
			pos: index,
			nucleotide,
			mutationPoss,
			proteinRegion,
		};
	});

	return [possibilityMap, windowSliceArr];
}

export const proteinRegionLabels = new Array(29904).map((_, idx) =>
	determineProteinRegion(idx)
);

export function determineProteinRegion(index, decimateFactor = 1) {
	for (const [protein, range] of Object.entries(proteinRegions)) {
		const [start, end] = range.split("-").map(Number);

		if (index * decimateFactor >= start && index * decimateFactor <= end) {
			return `${protein}`;
		}
	}
	return "";
}

const calculateMidPoints = (regions) => {
	return Object.entries(regions).map((item) => {
		const [name, range] = item;
		const [start, end] = range.split("-").map(Number);
		return {
			label: name,
			position: parseInt(start + (end - start) / 2),
		};
	});
};

export const proteinMidPoints = calculateMidPoints(proteinRegions);

function calculateMutationPoss(windowSlice) {
	const mutationPoss = { A: 0, T: 0, G: 0, C: 0 };

	for (const key in mutationPoss) {
		mutationPoss[key] = fooModelOutputGenerator();
	}

	return mutationPoss;
}

// input: mutationPoss and indexes
// output: Sum of the mutation possibilities for each protein region
// i.e: "S": 17,
//      "ORF1a" : 23
const calculateSumOfMutationProbability = (currentProteinRegion) => {
	return currentProteinRegion.reduce((acc, curr) => {
		const sum = Object.keys(curr.mutationPoss)
			.filter((key) => key !== curr.nucleotide)
			.reduce(
				(acc, key) => parseFloat(acc) + parseFloat(curr.mutationPoss[key]),
				0
			)
			.toFixed(2);
		return (parseFloat(acc) + parseFloat(sum)).toFixed(2);
	}, 0);
};

export function generateProteinRegionPossibility2(data) {
	const proteinRegionPossMap = {};
	const newData = [];
	for (const [proteinName, interval] of Object.entries(proteinRegions)) {
		const [start, end] = interval.split("-").map((pos) => parseInt(pos));
		const length = end - start;
		const currentProteinRegion = data.slice(start, end + 1);
		const prSum = parseFloat(
			calculateSumOfMutationProbability(currentProteinRegion)
		);
	}
}

export function generateProteinRegionPossibility(data) {
	// console.log("protein region");
	// console.log("data2: ", data);
	// console.log("map:", proteinRegions);
	const proteinRegionPossMap = {};
	for (const [proteinName, interval] of Object.entries(proteinRegions)) {
		const [start, end] = interval.split("-").map((pos) => parseInt(pos));
		const length = end - start;
		const currentProteinRegion = data.slice(start, end + 1);
		const prSum = parseFloat(
			calculateSumOfMutationProbability(currentProteinRegion)
		);

		proteinRegionPossMap[proteinName] = prSum;
	}
	console.log("proteinRegionPossMap:", proteinRegionPossMap);
	return proteinRegionPossMap;
}

// Example usage:
export function main(randomGenome) {
	const [possibilityMap, windowSliceArr] = windowPropagation(
		randomGenome,
		windowSize
	);
	// console.log("result:", possibilityMap);
	return [possibilityMap, windowSliceArr];
}
