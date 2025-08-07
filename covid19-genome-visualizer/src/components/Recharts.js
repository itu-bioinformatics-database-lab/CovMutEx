import React, { useEffect, useMemo, useRef, useState } from "react";
import { useSelector } from "react-redux";
import Chart from "chart.js/auto";
import annotationPlugin from "chartjs-plugin-annotation";
import zoomPlugin from "chartjs-plugin-zoom";
import { nucleotides } from "../helpers/helperFunctions";
import { proteinRegions } from "../data/proteinRegions";
import {
  proteinRegionColorMap,
  proteinRegionColorMapAnnotations,
} from "../utils/proteinRegionColorMap";
import SidePanel from "./SidePanel";

const nucleotideColors = {
  'A': 'rgba(240, 80, 80)',     // Red
  'C': 'rgba(80, 80, 240)',     // Blue
  'T': 'rgba(80, 240, 80)',     // Green
  'G': 'rgba(240, 200, 80)',   // Yellow
  'N': 'rgba(200, 160, 60, 1)'  // Gray (unknown)
};


const MINIMUM_BARS = 50;
const HIGH_DETAIL_THRESHOLD = 1000; // Preserved for high-res view transition

// Register plugins
Chart.register(
  annotationPlugin,
  zoomPlugin,
    {
    id: 'highResLoadingOverlay',
    afterDraw: (chart) => {
      if (chart.highResLoading) {
        const { ctx, chartArea: { top, left, width, height } } = chart;
        ctx.save();
        ctx.fillStyle = 'white';
        ctx.fillRect(left, top, width, height);
        ctx.font = '20px Arial';
        ctx.fillStyle = '#000000';
        ctx.textAlign = 'center';
        ctx.fillText('High-Resolution loading...', left + width / 2, top + height / 2);
        ctx.restore();
      }
    }
  },
  {
    id: 'weblogoOverlay',
    afterDraw: (chart) => {
      const { ctx, chartArea: { top, left, width, height } } = chart;
      if (chart.weblogoLoading) {
        ctx.save();
        ctx.fillStyle = '#FFFFFF';
        ctx.fillRect(left, top, width, height);
        ctx.font = '20px Arial';
        ctx.fillStyle = '#000000';
        ctx.textAlign = 'center';
        ctx.fillText('WebLogo Loading...', left + width / 2, top + height / 2);
        ctx.restore();
        return;
      }
      if (chart.weblogoImage && chart.weblogoMode && chart.weblogoImage.img) {
        const aspectRatio = chart.weblogoImage.width / chart.weblogoImage.height;
        let drawWidth = width;
        let drawHeight = width / aspectRatio;
        if (drawHeight > height) {
          drawHeight = height;
          drawWidth = height * aspectRatio;
        }
        const x = left + (width - drawWidth) / 2;
        const y = top + (height - drawHeight) / 2;
        ctx.drawImage(chart.weblogoImage.img, x, y, drawWidth, drawHeight);
      }
    }
  },
  {
    id: 'chartVisibilityController',
    beforeDraw: (chart) => {
      if (chart.weblogoImage && chart.weblogoMode) {
        const {ctx, chartArea: {top, left, width, height}} = chart;
        if (!chart.clearedForWeblogo) {
          ctx.save();
          ctx.globalCompositeOperation = 'destination-over';
          ctx.fillStyle = 'white';
          ctx.fillRect(left, top, width, height);
          ctx.restore();
          chart.clearedForWeblogo = true;
        }
      } else {
        chart.clearedForWeblogo = false;
      }
    }
  }
);


const GenomeChart = ({ genomeData, genomeSequence }) => {
  const chartRef = useRef(null);
  const viewRangeToPreserve = useRef(null);
  const [activeProtein, setActiveProtein] = useState(null);
  const [showFullAnnotation, setShowFullAnnotation] = useState(false);
  const [weblogoLoading, setWeblogoLoading] = useState(false);

  const selectedProteinRegion = useSelector(
    (state) => state.genome.selectedProteinRegion
  );
  
  const [focusedProtein, setFocusedProtein] = useState(null);

  useEffect(() => {
    setFocusedProtein(selectedProteinRegion);
  }, [selectedProteinRegion]);

  const [decimateFactor, setDecimateFactor] = useState(25);
  const [zoomLevel, setZoomLevel] = useState(30000);
  const [highResViewRange, setHighResViewRange] = useState(null);
  const currentDecimateFactor = highResViewRange ? 1 : decimateFactor;

  // *** FINAL FIX LOCATION ***
  // This useEffect correctly sets the initial view parameters when the focused protein changes.
  useEffect(() => {
    if (focusedProtein && proteinRegions[focusedProtein]) {
      const [start, end] = proteinRegions[focusedProtein].split("-").map(Number);
      const regionLength = end - start;
      
      // Ensure we display a minimum number of bars for readability.
      const newDecimateFactor = Math.max(1, Math.floor(regionLength / MINIMUM_BARS));
      
      setDecimateFactor(newDecimateFactor);
      setHighResViewRange(null);
      setZoomLevel(regionLength);
      viewRangeToPreserve.current = { min: start, max: end };
    } else if (!focusedProtein && !highResViewRange) {
      // Reset to default full-genome view
      setDecimateFactor(25);
      setHighResViewRange(null);
      setZoomLevel(30000);
      viewRangeToPreserve.current = null;
    }
  }, [focusedProtein]);
  
  const API_URL = process.env.REACT_APP_API_URL;
  
  const normalizedData = useMemo(() => {
    if (!genomeData || !genomeData.length) return [];
    const isAlreadyNormalized = Math.abs(nucleotides.reduce((sum, _, i) => sum + (genomeData[i]?.[0] || 0), 0) - 1.0) < 1e-9;
    if (isAlreadyNormalized) return genomeData;
    
    return nucleotides.map((_, nucIdx) =>
      genomeData[nucIdx].map((count, pos) => {
        const total = nucleotides.reduce((sum, _, i) => sum + (genomeData[i]?.[pos] || 0), 0);
        return total > 0 ? count / total : 0;
      })
    );
  }, [genomeData]);

  const chartViewData = useMemo(() => {
    if (focusedProtein && proteinRegions[focusedProtein] && normalizedData?.length) {
      const [start, end] = proteinRegions[focusedProtein].split('-').map(Number);
      const slicedData = normalizedData.map(arr => arr.slice(start, end + 1));
      return { dataForView: slicedData, offsetForView: start };
    }
    return { dataForView: normalizedData, offsetForView: 0 };
  }, [normalizedData, focusedProtein]);
  
  // All other functions (createAnnotations, getFullResolutionWebLogoData, handlers, etc.)
  // are correct as of the previous step. They are included here for completeness.
  const createAnnotations = () => {
    const chartInstance = chartRef.current?.chartInstance;
    if (chartInstance?.weblogoMode || (focusedProtein && !showFullAnnotation && !activeProtein)) { return []; }
    const mapPos = (genomePosition) => {
      if (highResViewRange) { return genomePosition - highResViewRange.min; }
      const startOffset = chartViewData.offsetForView;
      return (genomePosition - startOffset) / currentDecimateFactor;
    };
    if (activeProtein) {
      const [start, end] = proteinRegions[activeProtein].split("-").map(Number);
      return [{ display: true, type: "box", xMin: mapPos(start), xMax: mapPos(end), yMin: 0, yMax: 3, backgroundColor: proteinRegionColorMapAnnotations[activeProtein], borderColor: proteinRegionColorMap[activeProtein], borderWidth: 2, label: { content: activeProtein, enabled: true, position: "start" }, z: 10 }];
    } else if (showFullAnnotation) {
      return Object.keys(proteinRegions).map((key) => {
        const [start, end] = proteinRegions[key].split("-").map(Number);
        return { display: true, type: "box", xMin: mapPos(start), xMax: mapPos(end), yMin: 0, yMax: 3, backgroundColor: proteinRegionColorMapAnnotations[key], borderColor: proteinRegionColorMap[key], borderWidth: 2, label: { content: key, enabled: true, position: "start" }, z: 10 };
      });
    }
    return [];
  };
  
  const getFullResolutionWebLogoData = (startPos, endPos) => {
    if (!genomeData || !genomeData.length || !genomeSequence) return null;
    const numPositions = endPos - startPos + 1;
    const probabilityMatrix = Array.from({ length: numPositions }, () => [0, 0, 0, 0]);
    const referenceSequence = [];
    const positions = [];
    for (let i = 0; i < numPositions; i++) {
      const absolutePos = startPos + i;
      positions.push(absolutePos);
      referenceSequence.push(genomeSequence[absolutePos] || "N");
      const dataAccessIndex = absolutePos;
      let probabilityData = [0, 0, 0, 0];
      if (dataAccessIndex >= 0 && dataAccessIndex < (genomeData[0]?.length || 0)) {
          for (let nucIndex = 0; nucIndex < 4; nucIndex++) {
              probabilityData[nucIndex] = genomeData[nucIndex]?.[dataAccessIndex] || 0;
          }
      }
      const sum = probabilityData.reduce((a, b) => a + b, 0);
      if (sum > 0 && Math.abs(sum - 1.0) > 1e-9) {
          probabilityMatrix[i] = probabilityData.map(val => val / sum);
      } else {
          probabilityMatrix[i] = probabilityData;
      }
    }
    return { positions, probabilityMatrix, referenceSequence };
  };

  const fetchWebLogoImage = async (start, end) => {
    if (!genomeData || !genomeSequence) return null;
    const weblogoData = getFullResolutionWebLogoData(start, end);
    if (!weblogoData) return null;
    try {
      const response = await fetch(`${API_URL}/generate-weblogo/`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ start, end, probability_matrix: weblogoData.probabilityMatrix, reference_sequence: weblogoData.referenceSequence.join(''), nucleotide_order: nucleotides, decimate_factor: 1 })
      });
      if (!response.ok) { throw new Error(`HTTP ${response.status}`); }
      return URL.createObjectURL(await response.blob());
    } catch (error) { console.error("WebLogo generation failed:", error); return null; }
  };
  
  const displayWebLogo = async (chart, startPos, endPos) => {
    if (chart.weblogoTransition) return;
    if (focusedProtein) {
      const [proteinStart, proteinEnd] = proteinRegions[focusedProtein].split('-').map(Number);
      if (startPos < proteinStart || endPos > proteinEnd) { return; }
    }
    let controller;
    try {
      chart.weblogoTransition = true;
      controller = new AbortController();
      setWeblogoLoading(true);
      if (chart.weblogoImage) { URL.revokeObjectURL(chart.weblogoImage.url); chart.weblogoImage = null; }
      chart.weblogoLoading = true;
      chart.update();
      const imageUrl = await fetchWebLogoImage(startPos, endPos, controller.signal);
      if (!imageUrl) { chart.weblogoMode = false; chart.weblogoLoading = false; return; }
      const img = new Image();
      await new Promise((resolve, reject) => {
        img.onload = () => {
          chart.weblogoImage = { url: imageUrl, img, start: startPos, end: endPos, width: img.width, height: img.height };
          chart.weblogoLoading = false; chart.weblogoMode = true; chart.clearedForWeblogo = false; resolve();
        };
        img.onerror = () => reject(new Error('WebLogo image failed to load'));
        img.src = imageUrl;
      });
      chart.update();
    } catch (error) {
      if (error.name !== 'AbortError') {
        if (chart.weblogoImage) { URL.revokeObjectURL(chart.weblogoImage.url); chart.weblogoImage = null; }
        chart.weblogoMode = false; chart.weblogoLoading = false; chart.update();
      }
    } finally {
      if (controller) controller.abort();
      chart.weblogoTransition = false;
      setWeblogoLoading(false);
    }
  };

  const updateBarChart = (chart) => {
    if (chart.weblogoImage) {
      URL.revokeObjectURL(chart.weblogoImage.url);
      chart.weblogoImage = null;
      chart.weblogoMode = false;
      chart.update();
    }
  };

  const ZOOM_OUT_FIXED_BAR_COUNT = 60;
  const WEBLOGO_THRESHOLD = 450;
  const ZOOM_FACTOR = 1.30;

  const handleZoomIn = () => {
    const chart = chartRef.current?.chartInstance;
    if (!chart || chart.weblogoMode) return;
    let startOffset = highResViewRange ? highResViewRange.min : chartViewData.offsetForView;
    const effectiveDecimate = highResViewRange ? 1 : currentDecimateFactor;
    const { min: currentMin, max: currentMax } = chart.scales.x;
    const minGenomePos = startOffset + (currentMin * effectiveDecimate);
    const maxGenomePos = startOffset + (currentMax * effectiveDecimate);
    const currentRangeInBasePairs = maxGenomePos - minGenomePos;
    const newRangeInBasePairs = currentRangeInBasePairs / ZOOM_FACTOR;
    if (!highResViewRange && newRangeInBasePairs <= HIGH_DETAIL_THRESHOLD) {
      chart.highResLoading = true;
      chart.update('none');
      setTimeout(() => {
          const centerGenomePos = (minGenomePos + maxGenomePos) / 2;
          const newVisibleRange = 100;
          let newMinGenomePos = Math.round(centerGenomePos - newVisibleRange / 2);
          let newMaxGenomePos = newMinGenomePos + newVisibleRange;
          const [minBound, maxBound] = focusedProtein ? proteinRegions[focusedProtein].split('-').map(Number) : [0, 30000];
          newMinGenomePos = Math.max(minBound, newMinGenomePos);
          newMaxGenomePos = Math.min(maxBound, newMaxGenomePos);
          setHighResViewRange({ min: newMinGenomePos, max: newMaxGenomePos });
      }, 500);
      return;
    }
    if ((highResViewRange || decimateFactor === 1) && newRangeInBasePairs <= WEBLOGO_THRESHOLD) {
      const centerGenomePos = Math.round((minGenomePos + maxGenomePos) / 2);
      const [minBound, maxBound] = focusedProtein ? proteinRegions[focusedProtein].split('-').map(Number) : [1, 30000];
      let startPos = Math.max(minBound, centerGenomePos - 12);
      let endPos = Math.min(maxBound, startPos + 25);
      displayWebLogo(chart, startPos, endPos);
      return;
    }
    chart.zoom(ZOOM_FACTOR);
  };

  const handleZoomOut = () => {
    const chart = chartRef.current?.chartInstance;
    if (!chart) return;
    if (chart.weblogoMode) { updateBarChart(chart); return; }
    if (highResViewRange) {
      if (focusedProtein) {
        const [proteinStart, proteinEnd] = proteinRegions[focusedProtein].split('-').map(Number);
        const regionLength = proteinEnd - proteinStart;
        const newDecimateFactor = Math.max(1, Math.floor(regionLength / MINIMUM_BARS));
        viewRangeToPreserve.current = { min: proteinStart, max: proteinEnd };
        setHighResViewRange(null);
        setDecimateFactor(newDecimateFactor);
        setZoomLevel(regionLength);
      } else {
        const { min: currentMin, max: currentMax } = chart.scales.x;
        const minGenomePos = highResViewRange.min + currentMin;
        const maxGenomePos = highResViewRange.min + currentMax;
        const newGenomicRange = ZOOM_OUT_FIXED_BAR_COUNT * 25;
        const centerGenomePos = (minGenomePos + maxGenomePos) / 2;
        const newMinGenomePos = Math.max(0, Math.floor(centerGenomePos - newGenomicRange / 2));
        const newMaxGenomePos = newMinGenomePos + newGenomicRange;
        viewRangeToPreserve.current = { min: newMinGenomePos, max: newMaxGenomePos };
        setHighResViewRange(null);
        setDecimateFactor(25);
        setZoomLevel(newMaxGenomePos - newMinGenomePos);
      }
      return;
    }
    chart.zoom(1 / ZOOM_FACTOR);
  };

  const handlePan = (direction) => {
    const chart = chartRef.current?.chartInstance;
    if (!chart) return;
    if (chart.weblogoMode && chart.weblogoImage) {
      const panStep = 25;
      const [minBound, maxBound] = focusedProtein ? proteinRegions[focusedProtein].split('-').map(Number) : [1, 30000];
      let startPos = chart.weblogoImage.start + (direction * panStep);
      startPos = Math.max(minBound, Math.min(maxBound - 25, startPos));
      const endPos = Math.min(maxBound, startPos + 25);
      if (Math.abs(startPos - chart.weblogoImage.start) > 1) { displayWebLogo(chart, startPos, endPos); }
      return;
    }
    if (highResViewRange) {
        const currentRange = highResViewRange.max - highResViewRange.min;
        const panAmount = Math.round(currentRange * 0.25);
        let newMin = highResViewRange.min + (direction * panAmount);
        let newMax = highResViewRange.max + (direction * panAmount);
        const [minBound, maxBound] = focusedProtein ? proteinRegions[focusedProtein].split('-').map(Number) : [0, 30000];
        if (newMin < minBound) { newMin = minBound; newMax = minBound + currentRange; }
        if (newMax > maxBound) { newMax = maxBound; newMin = maxBound - currentRange; }
        setHighResViewRange({ min: Math.round(newMin), max: Math.round(newMax) });
        return;
    }
    const { min: currentMin, max: currentMax } = chart.scales.x;
    const viewRange = currentMax - currentMin;
    const panAmount = viewRange * direction * 0.25;
    const genomeLength = focusedProtein ? proteinRegions[focusedProtein].split('-').map(Number)[1] - proteinRegions[focusedProtein].split('-').map(Number)[0] : 30000;
    const genomeLengthChart = genomeLength / currentDecimateFactor;
    let newMin = currentMin + panAmount;
    let newMax = currentMax + panAmount;
    if (newMin < 0) { newMin = 0; newMax = viewRange; }
    else if (newMax > genomeLengthChart) { newMax = genomeLengthChart; newMin = newMax - viewRange; }
    chart.options.scales.x.min = newMin;
    chart.options.scales.x.max = newMax;
    chart.update('none');
  };

  const handleProteinRegionClick = (proteinName) => {
    if (focusedProtein === proteinName) { setFocusedProtein(null); } 
    else { setFocusedProtein(proteinName); }
  };

  const handleResetView = () => {
    const chart = chartRef.current?.chartInstance;
    if (chart?.weblogoMode) { updateBarChart(chart); }
    setFocusedProtein(null);
    setHighResViewRange(null);
  };

  useEffect(() => {
    const ctx = chartRef.current?.getContext("2d");
    const { dataForView, offsetForView } = chartViewData;

    if (!ctx || !dataForView?.length) return;
    if (chartRef.current.chartInstance) { chartRef.current.chartInstance.destroy(); }
    
    const chartDataSource = dataForView;
    const chartDataOffset = offsetForView;

    let dataForChart, dataStartPosition, dataEffectiveDecimateFactor;

    if (highResViewRange) {
        const { min: viewMin, max: viewMax } = highResViewRange;
        const relativeStart = viewMin - chartDataOffset;
        const relativeEnd = viewMax - chartDataOffset;
        const sliceStart = Math.max(0, relativeStart);
        const sliceEnd = Math.min(chartDataSource[0]?.length || 0, relativeEnd + 1);
        dataForChart = chartDataSource.map(dataset => dataset.slice(sliceStart, sliceEnd));
        dataStartPosition = chartDataOffset + sliceStart;
        dataEffectiveDecimateFactor = 1;
    } else {
        dataStartPosition = chartDataOffset;
        dataEffectiveDecimateFactor = decimateFactor;
        dataForChart = chartDataSource.map(dataset => {
            const decimated = [];
            for (let i = 0; i < dataset.length; i += dataEffectiveDecimateFactor) {
                const slice = dataset.slice(i, i + dataEffectiveDecimateFactor);
                if (slice.length > 0) { decimated.push(slice.reduce((a, b) => a + b, 0) / slice.length); }
            }
            return decimated;
        });
    }
    
    if (!dataForChart || !dataForChart[0]?.length) {
      return;
    }

    const decimatedLabels = Array.from({ length: dataForChart[0].length }, (_, idx) => {
        const position = dataStartPosition + (idx * dataEffectiveDecimateFactor);
        const nucleotide = genomeSequence[position] || "N";
        return `${position}-${nucleotide}`;
    });

    const mutationDetails = [];
    decimatedLabels.forEach((label, idx) => {
      const segmentStartPos = dataStartPosition + (idx * dataEffectiveDecimateFactor);
      const segmentEndPos = Math.min(30000, segmentStartPos + dataEffectiveDecimateFactor);
      const mutationProbs = { A: 0, T: 0, G: 0, C: 0 };
      let validPositions = 0;
      const refNucleotide = genomeSequence[segmentStartPos];
      const refIndex = nucleotides.indexOf(refNucleotide);
      if (refIndex === -1) return;
      for (let absolutePos = segmentStartPos; absolutePos < segmentEndPos; absolutePos++) {
        const dataAccessIndex = absolutePos;
        if (dataAccessIndex < 0 || dataAccessIndex >= (genomeData[0]?.length || 0)) continue;
        for (let targetIndex = 0; targetIndex < 4; targetIndex++) {
            if (targetIndex !== refIndex) {
                const targetNuc = nucleotides[targetIndex];
                mutationProbs[targetNuc] += genomeData[targetIndex]?.[dataAccessIndex] || 0;
            }
        }
        validPositions++;
      }
      if (validPositions > 0) { nucleotides.forEach(nuc => { mutationProbs[nuc] /= validPositions; }); }
      const totalMutProb = nucleotides.filter(nuc => nuc !== refNucleotide).reduce((sum, nuc) => sum + mutationProbs[nuc], 0);
      mutationDetails.push({ absolutePosition: segmentStartPos, refNuc: refNucleotide, mutations: { ...mutationProbs }, total: totalMutProb, });
    });
    
    const mutationDatasets = [];
    nucleotides.forEach(targetNuc => {
      const data = mutationDetails.map((detail) => (detail.refNuc === targetNuc || detail.total === 0) ? 0 : detail.mutations[targetNuc] || 0);
      mutationDatasets.push({ label: `${targetNuc} `, data, stack: 'mutation', maxBarThickness: 50, backgroundColor: (ctx) => ctx.chart.weblogoMode ? 'rgba(0,0,0,0)' : nucleotideColors[targetNuc], borderColor: (ctx) => ctx.chart.weblogoMode ? 'rgba(0,0,0,0)' : nucleotideColors[targetNuc] });
    });
    const filteredDatasets = mutationDatasets.filter(ds => ds.data.some(v => v > 0));
    const data = { labels: decimatedLabels, datasets: filteredDatasets };

    const options = {
      animation: false, responsive: true, maintainAspectRatio: false,
      scales: {
        x: { stacked: true, grid: { color: (c) => c.chart.weblogoMode ? 'rgba(0,0,0,0)' : 'rgba(0,0,0,0.1)' }, ticks: { color: (c) => c.chart.weblogoMode ? 'rgba(0,0,0,0)' : '#666' } },
        y: { stacked: true, type: 'logarithmic', min: 0.001, max: 1.0,title: {
    display: true,
    text: 'Mutation Probability (Log Scale)',
    color: (c) => c.chart.weblogoMode ? 'rgba(0,0,0,0)' : '#333'
  }, ticks: { callback: (v) => { if (v===1) return '10⁰'; if (v===0.1) return '10⁻¹'; if (v===0.01) return '10⁻²'; if (v===0.001) return '10⁻³'; return ''; }, color: (c) => c.chart.weblogoMode ? 'rgba(0,0,0,0)' : '#666' }, grid: { color: (c) => c.chart.weblogoMode ? 'rgba(0,0,0,0)' : 'rgba(0,0,0,0.1)' } },
      },
      plugins: {
        tooltip: {
          enabled: (c) => !c.chart.weblogoMode,
          callbacks: {
            title: (items) => `Position ${items[0].label.split("-")[0]} (Ref: ${items[0].label.split("-")[1]})`,
            label: () => null, 
            afterBody: (items) => {
              const detail = chartRef.current.chartInstance.mutationDetails?.[items[0].dataIndex];
              if (!detail) return [];
              const mutations = ['A','T','G','C'].filter(n => n !== detail.refNuc && detail.mutations[n] > 0.0001).map(n => `${detail.refNuc} → ${n}: ${detail.mutations[n].toFixed(4)}`);
              return [...mutations, `Total: ${detail.total.toFixed(4)}`];
            }
          }, displayColors: false
        },
        legend: { display: (c) => !c.chart.weblogoMode, },
        zoom: {
          zoom: { wheel: { enabled: true }, pinch: { enabled: true }, mode: "x",
            onZoomComplete: ({ chart }) => {
                let startOffset = highResViewRange ? highResViewRange.min : chartViewData.offsetForView;
                const effectiveDecimate = highResViewRange ? 1 : currentDecimateFactor;
                const { min: minIndex, max: maxIndex } = chart.scales.x;
                const minGenomePos = startOffset + (minIndex * effectiveDecimate);
                const maxGenomePos = startOffset + (maxIndex * effectiveDecimate);
                const newZoomLevel = maxGenomePos - minGenomePos;
                if (chart.weblogoMode && newZoomLevel > zoomLevel) { updateBarChart(chart); setZoomLevel(WEBLOGO_THRESHOLD * 1.2); return; }
                setZoomLevel(newZoomLevel);
                if (newZoomLevel <= HIGH_DETAIL_THRESHOLD && !highResViewRange) {
                    chart.highResLoading = true; chart.update('none');
                    setTimeout(() => {
                        const center = (minGenomePos + maxGenomePos) / 2; const range = 100;
                        const [minB, maxB] = focusedProtein ? proteinRegions[focusedProtein].split('-').map(Number) : [0, 30000];
                        let newMin = Math.max(minB, Math.round(center - range / 2));
                        let newMax = Math.min(maxB, newMin + range);
                        setHighResViewRange({ min: newMin, max: newMax });
                    }, 0);
                    return; 
                }
                if (newZoomLevel > HIGH_DETAIL_THRESHOLD && highResViewRange) {
                    const newRange = ZOOM_OUT_FIXED_BAR_COUNT * 25; const center = (minGenomePos + maxGenomePos) / 2;
                    let newMin = Math.floor(center - newRange / 2); let newMax = newMin + newRange;
                    viewRangeToPreserve.current = { min: Math.max(0, newMin), max: newMax };
                    setHighResViewRange(null); setDecimateFactor(25);
                    return;
                }
                if (newZoomLevel <= WEBLOGO_THRESHOLD && !chart.weblogoMode) {
                    const centerPos = (minGenomePos + maxGenomePos) / 2;
                    const [minB, maxB] = focusedProtein ? proteinRegions[focusedProtein].split('-').map(Number) : [1, 30000];
                    let startPos = Math.max(minB, Math.round(centerPos) - 12);
                    let endPos = Math.min(maxB, startPos + 25);
                    displayWebLogo(chart, startPos, endPos);
                }
            }
          },
          pan: { enabled: true, mode: "x" },
        },
        annotation: { annotations: createAnnotations(), },
      },
    };

    const chartInstance = new Chart(ctx, { type: "bar", data, options });
    chartInstance.mutationDetails = mutationDetails;
    chartRef.current.chartInstance = chartInstance;
    chartInstance.highResLoading = false;
    
    if (viewRangeToPreserve.current) {
        const { min, max } = viewRangeToPreserve.current;
        const startOffset = chartViewData.offsetForView;
        chartInstance.options.scales.x.min = (min - startOffset) / currentDecimateFactor;
        chartInstance.options.scales.x.max = (max - startOffset) / currentDecimateFactor;
        viewRangeToPreserve.current = null;
    } else if (!highResViewRange && !focusedProtein) {
        chartInstance.resetZoom();
    }
    
    return () => { chartInstance.destroy(); };
  }, [chartViewData, genomeData, genomeSequence, decimateFactor, highResViewRange]);

  useEffect(() => {
    const chartInstance = chartRef.current?.chartInstance;
    if (chartInstance) {
      chartInstance.options.plugins.annotation.annotations = createAnnotations();
      chartInstance.update();
    }
  }, [activeProtein, showFullAnnotation]);

  const handleProteinHover = (protein) => setActiveProtein(protein);
  const handleProteinLeave = () => setActiveProtein(null);
  const handleShowFullAnnotation = () => setShowFullAnnotation((prev) => !prev);

  return (
    <div className="overflow-x-hidden">
      <div className="chart-container w-full flex bg-[#f6f7f9] py-5">
        <SidePanel
          proteinRegions={proteinRegions}
          onProteinHover={handleProteinHover}
          onProteinLeave={handleProteinLeave}
          handleShowFullAnnotation={handleShowFullAnnotation}
          onProteinClick={handleProteinRegionClick}
        />
        <div className="w-full relative">
          {(focusedProtein || highResViewRange) && (
            <button
              onClick={handleResetView}
              className="absolute top-8 left-8 z-20 bg-white hover:bg-gray-100 text-gray-800 font-semibold py-1 px-3 border border-gray-400 rounded-lg shadow-md"
              aria-label="Reset View"
            >
              Reset View
            </button>
          )}

          <div className="absolute right-0 top-1/2 transform -translate-y-1/2 z-10">
            <div className="flex flex-col items-center bg-white rounded-lg shadow-md overflow-hidden">
              <button onClick={handleZoomIn} className="w-8 h-8 flex items-center justify-center bg-white hover:bg-gray-100 border-b border-gray-200" aria-label="Zoom in">
                <svg xmlns="http://www.w3.org/2000/svg" className="h-4 w-4 text-gray-700" viewBox="0 0 20 20" fill="currentColor"><path fillRule="evenodd" d="M10 5a1 1 0 011 1v3h3a1 1 0 110 2h-3v3a1 1 0 11-2 0v-3H6a1 1 0 110-2h3V6a1 1 0 011-1z" clipRule="evenodd" /></svg>
              </button>
              <button onClick={handleZoomOut} className="w-8 h-8 flex items-center justify-center bg-white hover:bg-gray-100" aria-label="Zoom out">
                <svg xmlns="http://www.w3.org/2000/svg" className="h-4 w-4 text-gray-700" viewBox="0 0 20 20" fill="currentColor"><path fillRule="evenodd" d="M5 10a1 1 0 011-1h8a1 1 0 110 2H6a1 1 0 01-1-1z" clipRule="evenodd" /></svg>
              </button>
            </div>
          </div>
          <div className="absolute -bottom-6 left-1/2 transform -translate-x-1/2 flex items-center space-x-4 z-10">
            <span className="text-sm text-gray-600">View earlier</span>
            <button onClick={() => handlePan(-1)} className="bg-white p-2 rounded-full shadow-md hover:bg-gray-100">
              <svg xmlns="http://www.w3.org/2000/svg" className="h-6 w-6" fill="none" viewBox="0 0 24 24" stroke="currentColor"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M15 19l-7-7 7-7" /></svg>
            </button>
            <button onClick={() => handlePan(1)} className="bg-white p-2 rounded-full shadow-md hover:bg-gray-100">
              <svg xmlns="http://www.w3.org/2000/svg" className="h-6 w-6" fill="none" viewBox="0 0 24 24" stroke="currentColor"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5l7 7-7 7" /></svg>
            </button>
            <span className="text-sm text-gray-600">View later</span>
          </div>
          <canvas className="w-full h-[90vh] max-h-screen bg-white mt-6 pl-4 pr-8 py-2 rounded-xl shadow-md" ref={chartRef} />
        </div>
      </div>
    </div>
  );
};

export default GenomeChart;
