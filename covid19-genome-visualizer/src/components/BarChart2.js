import React, { useEffect, useRef, useState } from "react";
import Chart from "chart.js/auto";
import zoomPlugin from "chartjs-plugin-zoom";
import { useDispatch, useSelector } from "react-redux";
import { resetChart } from "../features/genome/genomeSlice";
import { proteinRegionColorMap } from "../utils/proteinRegionColorMap";
import { GrPowerReset } from "react-icons/gr";
import { Button } from "@material-tailwind/react";
import ZoomSlider from "./ZoomSlider";

Chart.register(zoomPlugin);

let nucleotides = { 0: "A", 1: "C", 2: "T", 3: "G" };

function BarChart2({ data, seq }) {
    const chartRef = useRef();
    const dispatch = useDispatch();
    const chartTitle = useSelector((state) => state.genome.chartTitle);
    const [_chart, setChart] = useState(null);

    // --- GÜVENLİ DECIMATE ---
    const decimateData = (data, factor) => {
        // Eğer data yoksa veya boşsa, boş dataset yapısı dön
        if (!data || !Array.isArray(data) || data.length === 0) return [[], [], [], []];

        const aggregated = [[], [], [], []];

        for (let i = 0; i < data.length; i += factor) {
            const slice = data.slice(i, i + factor);
            
            const sums = [0, 0, 0, 0];
            let count = 0;

            slice.forEach(item => {
                // mutationPoss kontrolü
                if (item && item.mutationPoss) {
                    sums[0] += item.mutationPoss.A || 0;
                    sums[1] += item.mutationPoss.C || 0;
                    sums[2] += item.mutationPoss.T || 0;
                    sums[3] += item.mutationPoss.G || 0;
                    count++;
                }
            });

            if (count > 0) { 
                 aggregated[0].push(sums[0] / count);
                 aggregated[1].push(sums[1] / count);
                 aggregated[2].push(sums[2] / count);
                 aggregated[3].push(sums[3] / count);
            } else {
                // Veri yoksa 0 bas
                aggregated[0].push(0);
                aggregated[1].push(0);
                aggregated[2].push(0);
                aggregated[3].push(0);
            }
        }
        return aggregated;
    };

    useEffect(() => {
        if (!data || !seq) return;

        const currentChartRef = chartRef.current;
        const decimatedData = decimateData(data, 15); 
        
        // Etiketler (seq uzunluğu kontrol edilerek)
        const labelsLength = decimatedData[0].length;
        const labels = seq.slice(0, labelsLength).split(""); 
        const label_indexes = labels.map((_, idx) => idx);

        const datasets = decimatedData.map((dataset, idx) => ({
            label: nucleotides[idx],
            data: dataset,
            borderColor: getColorForNucleotide(nucleotides[idx]),
            backgroundColor: getColorForNucleotide(nucleotides[idx]),
        }));

        if (currentChartRef) {
            if (_chart) _chart.destroy();

            const ctx = currentChartRef.getContext("2d");

            if (ctx) {
                const chart = new Chart(ctx, {
                    type: "bar",
                    data: {
                        labels: label_indexes,
                        datasets: datasets,
                    },
                    options: {
                        animation: false,
                        scales: { x: { stacked: true }, y: { stacked: true } },
                        responsive: true,
                        maintainAspectRatio: false,
                        plugins: {
                            zoom: {
                                zoom: {
                                    wheel: { enabled: true },
                                    pinch: { enabled: true },
                                    mode: "x",
                                },
                                pan: { enabled: true, mode: "x" },
                            },
                            title: {
                                display: true,
                                color: proteinRegionColorMap[chartTitle] || "#000000",
                                position: "bottom",
                                text: chartTitle || "Genome Mutation Risk",
                                font: { size: 16 },
                            },
                        },
                    },
                });
                setChart(chart);
                return () => chart.destroy();
            }
        }
    // eslint-disable-next-line
    }, [data, seq, chartTitle]);

    const handleReset = () => {
        dispatch(resetChart());
        if (_chart) _chart.resetZoom();
    };

    const handleZoom = (zoomLevel) => {
        if (_chart) _chart.zoom(zoomLevel);
    };

    return (
        <div className="flex">
            <div className="flex items-center justify-center chart-container">
                <ZoomSlider handleZoom={handleZoom} />
            </div>
            <div className="chart-container" style={{ height: "50vh", width: "60vw" }}>
                <div className="ml-12 mb-2">
                    <Button color="blue" variant="gradient" className="flex gap-2 justify-center items-center" onClick={handleReset}>
                        <div>Reset Chart</div>
                        <GrPowerReset size={20} color="white" />
                    </Button>
                </div>
                <canvas ref={chartRef} />
            </div>
        </div>
    );
}

function getColorForNucleotide(nucleotide) {
    const colorMap = { A: "#FF5733", C: "#3399FF", T: "#33CC33", G: "#9966FF" };
    return colorMap[nucleotide] || "#000000";
}

export default BarChart2;