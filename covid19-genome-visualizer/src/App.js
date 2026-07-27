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
import PayloadRenderer from "./components/PayloadRenderer";
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
  // Visualization is delegated to <PayloadRenderer/>, which reads the v2.0
  // predictionPayload (and legacy fields) directly from the store. We only
  // need a tiny slice here for the page-level guards.
  const {
    dataset: genomeData,              // [{mutationPoss}] format for BarChart
    isSelected,
    loading,
  } = useSelector((state) => state.genome);

  // Linear / logarithmic toggle applies only to the categorical ATGC chart.
  // PayloadRenderer forwards it to GenomeChart; non-categorical tracks ignore.
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
  const handleNavbarSubmit = async (inputParams) => {
    const params = {
      nodeId: inputParams?.nodeId || "default_node_id",
      elapsedDay: inputParams?.elapsedDay ? Number(inputParams.elapsedDay) : 0,
      selectedModel: inputParams?.selectedModel || "balanced_data_model",
      selectedProteinRegion: inputParams?.selectedProteinRegion || null,
      isNewUpload: false,
      customParameters: inputParams?.customParameters || {},
      variant: inputParams?.variant || undefined,
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

        {/* VISUAL COMPARE SAYFASI */}
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

              {/* CHARTS — dispatched by task.kind in the v2.0 PredictionPayload */}
              <PayloadRenderer scaleType={scaleType} />
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
