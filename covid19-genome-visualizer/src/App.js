import { useEffect, useState } from "react";
import { MdTimeline, MdShowChart } from "react-icons/md";
import { useDispatch, useSelector } from "react-redux";
import { Routes, Route, useLocation, useNavigate } from "react-router-dom";
import "./App.css";

// Redux Actions & Slices
import {
  fetchPrediction,
  resetProteinRegion,
  loadNodesAndModels,
} from "./features/genome/genomeSlice";

// Components
import Navbar from "./components/Navbar";
import GenomeChart from "./components/Recharts";
import DoughnutChart from "./components/DoughnutChart";
import Error from "./components/Error";
import Contact from "./components/Contact";
import Nav from "./components/Nav";
import { About } from "./components/About";
import UploadModel from "./components/UploadModel";
import BenchmarkDashboard from "./components/BenchmarkDashboard";
import CompareModels from "./components/CompareModels";

// Assets
import logo from "./CovMutexLogo-removebg-preview.png";

function App() {
  const dispatch = useDispatch();
  const navigate = useNavigate();
  const location = useLocation();

  // --- REDUX SELECTORS ---
  const {
    genomeDataRaw,                    // [4][N] format for GenomeChart
    dataset: genomeData,              // [{mutationPoss}] format for BarChart
    genome: genomeSequence,
    proteinMutationProbs: protein_mutation_probs,
    selectedProteinRegion,
    isSelected,
    loading,
  } = useSelector((state) => state.genome);

  const [scaleType, setScaleType] = useState("logarithmic");

  // --- INITIAL DATA LOADING ---
  useEffect(() => {
    dispatch(loadNodesAndModels());
  }, [dispatch]);

  // --- ROUTING & LIFECYCLE ---
  useEffect(() => {
    if (location.pathname === "/") {
      dispatch(resetProteinRegion());
    }
  }, [location.pathname, dispatch]);

  // Sayfa yenilemeyi önleme (Veri varken)
  useEffect(() => {
    const preventRefresh = (e) => {
      if (isSelected || (genomeData && genomeData.length > 0)) {
        e.preventDefault();
        e.returnValue = "";
      }
    };
    window.addEventListener("beforeunload", preventRefresh);
    return () => window.removeEventListener("beforeunload", preventRefresh);
  }, [isSelected, genomeData]);

  // --- HANDLERS ---
  const handleNavbarSubmit = async (
    nodeId,
    elapsedDay,
    selectedModel,
    selectedProteinRegion
  ) => {
    const params = {
      nodeId: nodeId || "default_node_id",
      elapsedDay: elapsedDay ? Number(elapsedDay) : 0,
      selectedModel: selectedModel || "balanced_data_model",
      selectedProteinRegion: selectedProteinRegion || null,
      isNewUpload: false,
    };

    try {
      await dispatch(fetchPrediction(params)).unwrap();
      navigate("/genome-mutation-visualization");
    } catch (error) {
      console.error("Prediction failed:", error);
      navigate("/error", { replace: true });
    }
  };

  return (
    <div className="overflow-y-hidden min-h-screen flex flex-col">
      <Nav />
      
      <Routes>
        {/* ANA SAYFA (Input Formu) */}
        <Route
          exact
          path="/"
          element={
            <Navbar 
              onNodeSelect={() => {}} 
              onSubmit={handleNavbarSubmit} 
              isLoading={loading}
            />
          }
        />

        {/* HATA SAYFASI */}
        <Route path="/error" element={<Error />} />

        {/* YENİ UPLOAD SAYFASI */}
        <Route path="/upload-model" element={<UploadModel />} />

        {/* BENCHMARK / COMPARISON SAYFASI */}
        <Route path="/benchmark" element={<BenchmarkDashboard />} />

        {/* VISUALIZATION SAYFASI */}
        <Route path="/compare" element={<CompareModels />} />

        {/* SONUÇ GÖRSELLEŞTİRME SAYFASI */}
        <Route
          exact
          path="/genome-mutation-visualization"
          element={
            <div className="bg-[#f6f7f9] relative">
              {/* HEADER */}
              <div className="flex items-center justify-between px-8 pt-4 pb-0">
                <div className="flex items-center gap-3">
                  <img
                    src={logo}
                    className="w-[5rem] h-auto"
                    alt="CovMutEx Logo"
                  />
                  <h1 className="font-bold text-xl text-gray-800">
                    Genome Sequence Mutation Visualization
                  </h1>
                </div>
                {/* Scale Toggle */}
                <div className="flex items-center gap-1 bg-white rounded-lg shadow-sm border border-gray-200 p-0.5">
                  <button
                    onClick={() => setScaleType("logarithmic")}
                    className={`flex items-center gap-1.5 px-3 py-1.5 rounded-md text-sm font-medium transition-colors ${
                      scaleType === "logarithmic"
                        ? "bg-blue-600 text-white shadow-sm"
                        : "text-gray-600 hover:bg-gray-100"
                    }`}
                  >
                    <MdTimeline size={16} />
                    Log Scale
                  </button>
                  <button
                    onClick={() => setScaleType("linear")}
                    className={`flex items-center gap-1.5 px-3 py-1.5 rounded-md text-sm font-medium transition-colors ${
                      scaleType === "linear"
                        ? "bg-blue-600 text-white shadow-sm"
                        : "text-gray-600 hover:bg-gray-100"
                    }`}
                  >
                    <MdShowChart size={16} />
                    Linear Scale
                  </button>
                </div>
              </div>

              {/* LOADING INDICATOR */}
              {loading && (
                <div className="absolute inset-0 bg-white/80 z-50 flex items-center justify-center">
                  <div className="text-xl font-semibold text-blue-600 animate-pulse">
                    Calculating Predictions...
                  </div>
                </div>
              )}

              {/* CHARTS */}
              <div className="block">
                {genomeDataRaw && genomeDataRaw.length > 0 && (
                  <GenomeChart
                    genomeData={genomeDataRaw}
                    genomeSequence={genomeSequence}
                    scaleType={scaleType}
                  />
                )}
              </div>
              {!selectedProteinRegion && protein_mutation_probs && (
                <div className="flex justify-center py-6">
                  <DoughnutChart data={protein_mutation_probs} />
                </div>
              )}
            </div>
          }
        />

        {/* DİĞER SAYFALAR */}
        <Route exact path="/contact-us" element={<Contact />} />
        <Route exact path="/about" element={<About />} />
      </Routes>
    </div>
  );
}

export default App;