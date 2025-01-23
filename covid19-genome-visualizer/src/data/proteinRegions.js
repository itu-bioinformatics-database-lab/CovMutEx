export const proteinRegions = {
	ORF1ab: "266-21555",
	S: "21563-25384",
	ORF3a: "25393-26220",
	E: "26245-26472",
	M: "26523-27191",
	ORF6: "27202-27387",
	ORF7a: "27394-27759",
	ORF7b: "27756-27887",
	ORF8: "27894-28259",
	N: "28274-29533",
	ORF10: "29558-29674",
};

// export const proteinRegions = {
// 	ORF1ab: "0-399",
// 	S: "400-459",
// 	ORF3a: "460-510",
// 	E: "511-600",
// 	M: "601-630",
// 	ORF6: "631-699",
// 	ORF7a: "700-799",
// 	ORF8: "800-849",
// 	N: "850-899",
// 	ORF10: "900-999",
// };

export const proteinRegionsSize = getProteinRegionSize();

function getProteinRegionSize() {
	let a = {};
	Object.entries(proteinRegions).map((item) => {
		const [protein, range] = item;
		const [start, end] = range.split("-").map(Number);
		const length = end - start + 1;
		a[protein] = length;
		return {
			protein,
			length,
		};
	});
	return a;
}

// ORF1ab: { interval: "0-399", color: "rgba(255, 99, 132, 0.5)" },
// S: { interval: "400-459", color: "rgba(54, 102, 235, 0.5)" },
// ORF3a: { interval: "460-510", color: "rgba(255, 206, 86, 0.5)" },
// E: { interval: "511-600", color: "rgba(75, 192, 192, 0.5)" },
// M: { interval: "601-630", color: "rgba(153, 102, 255, 0.5)" },
// ORF6: { interval: "631-669", color: "rgba(255, 159, 64, 0.5)" },
// ORF7a: { interval: "700-799", color: "rgba(255, 99, 132, 0.5)" },
// ORF8: { interval: "800-849", color: "rgba(54, 162, 235, 0.5)" },
// N: { interval: "850-899", color: "rgba(255, 206, 86, 0.5)" },
// ORF10: { interval: "900-999", color: "rgba(75, 192, 192, 0.5)" },
