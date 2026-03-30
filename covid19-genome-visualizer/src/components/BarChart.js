import React, { useEffect, useRef, useState } from "react";
import Chart from "chart.js/auto";
import zoomPlugin from "chartjs-plugin-zoom";
import { useDispatch, useSelector } from "react-redux";
import { resetProteinRegion } from "../features/genome/genomeSlice"; 
import { proteinRegionColorMap } from "../utils/proteinRegionColorMap";

Chart.register(zoomPlugin);

function BarChart({ data }) {
    const chartRef = useRef();
    const dispatch = useDispatch();
    const [chunkView, setChunkView] = useState(true);
    const chartTitle = useSelector((state) => state.genome.chartTitle);
    const [_chart, setChart] = useState(null);

    // --- GÜVENLİ DECIMATE FONKSİYONU ---
    const decimateData = (data, factor) => {
        // Veri yoksa veya dizi değilse boş dizi dön
        if (!data || !Array.isArray(data) || data.length === 0) return [];
        
        return data.reduce((acc, _, index) => {
            if (index % factor === 0) {
                const slice = data.slice(index, index + factor);
                const averaged = slice.reduce((avg, curr) => {
                    // KONTROL: curr.mutationPoss var mı? Yoksa işlem yapma.
                    if (curr && curr.mutationPoss) {
                        Object.keys(curr.mutationPoss).forEach((key) => {
                            avg[key] = (avg[key] || 0) + curr.mutationPoss[key] / factor;
                        });
                    }
                    return avg;
                }, {});
                
                // Eğer averaged içinde veri oluştuysa listeye ekle
                if (Object.keys(averaged).length > 0) {
                    acc.push({
                        pos: index / factor,
                        nucleotide: "",
                        mutationPoss: averaged,
                    });
                }
            }
            return acc;
        }, []);
    };

    useEffect(() => {
        const currentChartRef = chartRef.current;
        const decimatedData = decimateData(data, 30); 

        // Eğer çizilecek veri yoksa dur
        if (decimatedData.length === 0) return;

        // İlk elemanı kontrol et
        const firstEntry = decimatedData[0];
        if (!firstEntry || !firstEntry.mutationPoss) return;

        const datasets = Object.keys(firstEntry.mutationPoss).map(
            (nucleotide) => ({
                label: nucleotide,
                // mutationPoss yoksa 0 ata
                data: decimatedData.map((entry) => entry.mutationPoss ? entry.mutationPoss[nucleotide] : 0),
                borderColor: getColorForNucleotide(nucleotide),
                backgroundColor: getColorForNucleotide(nucleotide),
            })
        );

        if (currentChartRef) {
            if (_chart) _chart.destroy();

            const ctx = currentChartRef.getContext("2d");

            if (ctx) {
                const labels = decimatedData.map(
                    (entry) => `${Math.floor(entry.pos)}-${entry.nucleotide}`
                );
               
                const chart = new Chart(ctx, {
                    type: "bar",
                    data: {
                        labels: labels,
                        datasets: datasets,
                    },
                    options: {
                        animation: false,
                        scales: {
                            x: {
                                type: "category",
                                stacked: true,
                            },
                            y: {
                                beginAtZero: true,
                                stacked: true,
                                max: 4,
                            },
                        },
                        plugins: {
                            zoom: {
                                zoom: {
                                    wheel: { enabled: true },
                                    pinch: { enabled: true },
                                    mode: "x",
                                },
                                pan: {
                                    enabled: true,
                                    mode: "x",
                                },
                            },
                            title: {
                                display: true,
                                color: proteinRegionColorMap[chartTitle] || "#000000",
                                position: "bottom",
                                text: chartTitle,
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
    }, [data, chartTitle]); 

    const handleReset = () => {
        setChunkView(true);
        dispatch(resetProteinRegion()); 
    };

    return (
        <div>
            <div className="flex justify-between items-center mb-2">
                <button onClick={handleReset} className="bg-gray-200 px-3 py-1 rounded text-sm hover:bg-gray-300">Reset View</button>
                <div className="flex gap-2">
                     <button className="bg-gray-200 px-2 py-1 rounded hover:bg-gray-300" onClick={() => _chart && _chart.zoom(1.1)}>+</button>
                     <button className="bg-gray-200 px-2 py-1 rounded hover:bg-gray-300" onClick={() => _chart && _chart.zoom(0.9)}>-</button>
                </div>
            </div>
            <canvas ref={chartRef} />
        </div>
    );
}

function getColorForNucleotide(nucleotide) {
    const colorMap = { A: "#FF5733", C: "#3399FF", T: "#33CC33", G: "#9966FF" };
    return colorMap[nucleotide] || "#000000";
}

export default BarChart;