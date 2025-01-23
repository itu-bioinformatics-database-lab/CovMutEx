export function getColorForNucleotide(nucleotide) {
  const colorMap = {
    A: "#FF5733",
    C: "#3399FF",
    T: "#33CC33",
    G: "#9966FF",
  };

  return colorMap[nucleotide] || "rgba(128, 128, 128, 0.6)";
}

export const nucleotides = ["A", "C", "T", "G"];

