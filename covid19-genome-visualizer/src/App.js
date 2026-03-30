import { useEffect, useState } from "react";
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
import ContextOverlay from "./components/ContextOverlay";

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

  // Context overlay params
  const [contextNodeId, setContextNodeId] = useState("");
  const [contextElapsedDay, setContextElapsedDay] = useState(60);
  const [contextProteinRegion, setContextProteinRegion] = useState("");

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

    // Store for context overlay
    setContextNodeId(nodeId || "");
    setContextElapsedDay(elapsedDay ? Number(elapsedDay) : 60);
    setContextProteinRegion(selectedProteinRegion || "");

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

        {/* VISUAL COMPARE SAYFASI */}
        <Route path="/compare" element={<CompareModels />} />

        {/* SONUÇ GÖRSELLEŞTİRME SAYFASI */}
        <Route
          exact
          path="/genome-mutation-visualization"
          element={
            <div className="bg-[#f6f7f9] relative min-h-screen">
              {/* HEADER / LOGO */}
              <h1 className="text-center pt-4 pb-0 font-bold text-xl text-gray-800">
                Genome Sequence Mutation Visualization
              </h1>
              <div className="absolute top-0 flex justify-center items-center">
                <img
                  src={logo}
                  className="w-[7rem] h-auto ml-[5.5rem]"
                  alt="CovMutEx Logo"
                />
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
              <div className="block md:flex md:justify-normal p-4">
                {/* GenomeChart uses genomeDataRaw which is [4][N] format */}
                {genomeDataRaw && genomeDataRaw.length > 0 && (
                  <div className="flex-1" style={{ height: "85vh", maxHeight: "85vh", overflow: "hidden" }}>
                    <GenomeChart
                      genomeData={genomeDataRaw}
                      genomeSequence={genomeSequence}
                    />
                  </div>
                )}

                {!selectedProteinRegion && protein_mutation_probs && (
                  <div className="md:w-1/3 mt-8 md:mt-0 flex justify-center">
                    <DoughnutChart data={protein_mutation_probs} />
                  </div>
                )}
              </div>

              {/* FR-4: Context Overlay */}
              <ContextOverlay
                nodeId={contextNodeId}
                elapsedDay={contextElapsedDay}
                selectedProteinRegion={contextProteinRegion}
              />
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