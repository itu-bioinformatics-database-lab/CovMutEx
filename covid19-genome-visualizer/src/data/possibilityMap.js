export function createArrayOfDictionaries() {
	const arraySize = 30000;
	return Array.from({ length: arraySize }, (_, index) => ({
		position: index,
		nucleotide: "",
		mutationPossibilities: {},
	}));
}

// Example usage:
