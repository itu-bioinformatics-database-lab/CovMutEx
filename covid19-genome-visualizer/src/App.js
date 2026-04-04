import { useEffect, useState } from "react";
import { useDispatch, useSelector } from "react-redux";
import { Routes, Route, useLocation, useNavigate } from "react-router-dom";
import { AnimatePresence, motion } from "framer-motion";
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

// Light page transition - fast fade only, no blocking
const PageTransition = ({ children }) => (
  <motion.div
    initial={{ opacity: 0 }}
    animate={{ opacity: 1 }}
    transition={{ duration: 0.15, ease: "easeOut" }}
  >
    {children}
  </motion.div>
);

function App() {
  const dispatch = useDispatch();
  const navigate = useNavigate();
  const location = useLocation();

  // --- REDUX SELECTORS ---
  const {
    genomeDataRaw,
    dataset: genomeData,
    genome: genomeSequence,
    proteinMutationProbs: protein_mutation_probs,
    selectedProteinRegion,
    isSelected,
    loading,
  } = useSelector((state) => state.genome);

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

  // Sayfa yenilemeyi onleme (Veri varken)
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
    <div className="min-h-screen flex flex-col bg-gray-50 dark:bg-gray-950 transition-colors duration-300">
      <Nav />

      <main className="flex-1">
        <AnimatePresence>
          <Routes location={location} key={location.pathname}>
            {/* ANA SAYFA */}
            <Route
              exact
              path="/"
              element={
                <PageTransition>
                  <Navbar
                    onNodeSelect={() => {}}
                    onSubmit={handleNavbarSubmit}
                    isLoading={loading}
                  />
                </PageTransition>
              }
            />

            {/* HATA SAYFASI */}
            <Route
              path="/error"
              element={
                <PageTransition>
                  <Error />
                </PageTransition>
              }
            />

            {/* UPLOAD SAYFASI */}
            <Route
              path="/upload-model"
              element={
                <PageTransition>
                  <UploadModel />
                </PageTransition>
              }
            />

            {/* BENCHMARK */}
            <Route
              path="/benchmark"
              element={
                <PageTransition>
                  <BenchmarkDashboard />
                </PageTransition>
              }
            />

            {/* VISUAL COMPARE */}
            <Route
              path="/compare"
              element={
                <PageTransition>
                  <CompareModels />
                </PageTransition>
              }
            />

            {/* SONUC GORSELLESTIRME - no animation wrapper for heavy chart */}
            <Route
              exact
              path="/genome-mutation-visualization"
              element={
                  <div className="bg-[#f6f7f9] dark:bg-gray-900 relative min-h-screen transition-colors">
                    {/* HEADER / LOGO */}
                    <h1 className="text-center pt-4 pb-0 font-bold text-xl text-gray-800 dark:text-gray-200">
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
                      <div className="absolute inset-0 bg-white/80 dark:bg-gray-900/80 z-50 flex items-center justify-center backdrop-blur-sm">
                        <div className="text-xl font-semibold text-blue-600 dark:text-blue-400 animate-pulse">
                          Calculating Predictions...
                        </div>
                      </div>
                    )}

                    {/* CHARTS */}
                    <div className="block md:flex md:justify-normal">
                      {genomeDataRaw && genomeDataRaw.length > 0 && (
                        <GenomeChart
                          genomeData={genomeDataRaw}
                          genomeSequence={genomeSequence}
                        />
                      )}

                      {!selectedProteinRegion && protein_mutation_probs && (
                        <DoughnutChart data={protein_mutation_probs} />
                      )}
                    </div>
                  </div>
              }
            />

            {/* DIGER SAYFALAR */}
            <Route
              exact
              path="/contact-us"
              element={
                <PageTransition>
                  <Contact />
                </PageTransition>
              }
            />
            <Route
              exact
              path="/about"
              element={
                <PageTransition>
                  <About />
                </PageTransition>
              }
            />
          </Routes>
        </AnimatePresence>
      </main>
    </div>
  );
}

export default App;
