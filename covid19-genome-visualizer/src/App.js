import { useDispatch, useSelector } from "react-redux";
import "./App.css";
import {
  loadNodesAndModels,
  setDataset,
  resetProteinRegion,
} from "./features/genome/genomeSlice";
import { useEffect, useState } from "react";
import Navbar from "./components/Navbar";
import GenomeChart from "./components/Recharts";
import DoughnutChart from "./components/DoughnutChart";
import ZoomChart from "./components/Test";
import { Routes, Route, useLocation, useNavigate } from "react-router-dom";
import Error from "./components/Error";
import logo from "./CovMutexLogo-removebg-preview.png";
import Contact from "./components/Contact";
import Nav from "./components/Nav";
import { About } from "./components/About";
import KnownHotspotCaseStudyPage from "./components/KnownHotspotCaseStudyPage";
import SpikeMutationPanel from "./components/SpikeMutationPanel";
 


function App() {
  const dispatch = useDispatch();
  const navigate = useNavigate();
  const location = useLocation();

  const nodeIds = useSelector((state) => state.genome.nodeList);
  const elapsedDay = useSelector((state) => state.genome.elapsedDay);
  const selectedModel = useSelector((state) => state.genome.model);
  const selectedProteinRegion = useSelector(
    (state) => state.genome.selectedProteinRegion
  );
  const isSelected = useSelector((state) => state.genome.isSelected);

  const [loading, setLoading] = useState(true);
  const [dataLoading, setDataLoading] = useState(false);
  const [genomeSequence, setGenomeSequence] = useState("");
  const [genomeData, setGenomeData] = useState([]);
  const [proteinRegionPossibilities, setProteinRegionPossibilities] = useState(
    {}
  );
  const [protein_mutation_probs, setProteinMutationProbs] = useState({});
  const [priestAnnotation, setPriestAnnotation] = useState(null);
  const [selectedAlgorithm, setSelectedAlgorithm] = useState(null);
  const API_URL = process.env.REACT_APP_API_URL;

  // Reset `selectedProteinRegion` on navigation
  useEffect(() => {
    if (location.pathname === "/") {
      dispatch(resetProteinRegion());
      setSelectedAlgorithm(null);
    }
  }, [location.pathname, dispatch]);

  // Prevent page refresh
  useEffect(() => {
    const preventRefresh = (e) => {
      if (isSelected || genomeData.length > 0) {
        e.preventDefault();
        e.returnValue = "";
      }
    };

    window.addEventListener("beforeunload", preventRefresh);

    return () => {
      window.removeEventListener("beforeunload", preventRefresh);
    };
  }, [isSelected, genomeData]);
//`${API_URL}/api/predict/`
// In App.js
const onSubmit = async (nodeId, elapsedDay, selectedModel, event) => {
  setDataLoading(true);
  setPriestAnnotation(null);
  setSelectedAlgorithm(selectedModel);
  try {
    const priestOnly = selectedModel === "PRIEST";
    // === FETCH 1: Get data based on user selection (for specific calcs) ===
    const response = await fetch(`${API_URL}/api/predict/`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        nodeId: nodeId || "default_node_id",
        elapsedDay: elapsedDay ? Number(elapsedDay) : 0,
        selectedModel: selectedModel || "default_model_path",
        selectedProteinRegion: selectedProteinRegion || null,
      }),
    });

    if (!response.ok) {
      navigate("/error", { replace: true });
      return;
    }

    const data = await response.json();

    // Store protein-specific calculations from the first fetch
    setProteinRegionPossibilities(data.proteinRegionPossibilities || {});
    setProteinMutationProbs(data.protein_mutation_probs || {});
    setPriestAnnotation(
      priestOnly
        ? {
            selected_node: data.selected_node,
            selected_node_name: data.selected_node_name,
            selected_node_accession: data.selected_node_accession,
            node_date: data.node_date,
            priest_period: data.priest_period,
            priest_period_available: data.priest_period_available,
            priest_score_method: data.priest_score_method,
            spike_mutations: data.spike_mutations || [],
            synonymous_spike_mutations: data.synonymous_spike_mutations || [],
            priest_summary: data.priest_summary || null,
          }
        : null
    );

    let finalGenomeData = data.genomeData;
    let finalGenomeSequence = data.genomeSequence;

    if (priestOnly) {
      setGenomeData([]);
      setGenomeSequence(data.genomeSequence || "");
      setProteinMutationProbs({});
      dispatch(
        setDataset({
          dataset: [],
          genome: data.genomeSequence || "",
          pr_poss: data.proteinRegionPossibilities || {},
          isSelected: true,
          protein_mutation_probs: {},
          selectedProteinRegion: null,
        })
      );
      navigate("/genome-mutation-visualization");
      return;
    }

    // === FETCH 2 (Conditional): If a region was selected, get the FULL genome ===
    // This ensures our chart always has the complete, correct data source.
    if (selectedProteinRegion) {
      console.log("Region selected. Fetching full genome data for the chart...");
      const fullGenomeResponse = await fetch(`${API_URL}/api/predict/`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          nodeId: nodeId || "default_node_id",
          elapsedDay: elapsedDay ? Number(elapsedDay) : 0,
          selectedModel: selectedModel || "default_model_path",
          selectedProteinRegion: null, // Ask for the full genome
        }),
      });
      const fullGenomeDataPayload = await fullGenomeResponse.json();
      finalGenomeData = fullGenomeDataPayload.genomeData;
      finalGenomeSequence = fullGenomeDataPayload.genomeSequence;
      console.log("Full genome data received, length:", finalGenomeData[0].length);
    }

    // Now, save the GUARANTEED full data to state and Redux
    setGenomeData(finalGenomeData);
    setGenomeSequence(finalGenomeSequence);

    dispatch(
      setDataset({
        dataset: finalGenomeData, // Pass the full dataset
        genome: finalGenomeSequence, // Pass the full sequence
        pr_poss: data.proteinRegionPossibilities,
        isSelected: true, // Mark as selected
        protein_mutation_probs: data.protein_mutation_probs,
        selectedProteinRegion: selectedProteinRegion || null,
      })
    );

    navigate("/genome-mutation-visualization");
  } catch (error) {
    console.error("Error during prediction:", error);
    navigate("/error", { replace: true });
  } finally {
    setDataLoading(false);
  }
};

  return (
    <div className="overflow-x-hidden">
      <Nav />
      <Routes>
        <Route
          exact
          path="/"
          element={<Navbar onNodeSelect={() => {}} onSubmit={onSubmit} />}
        />
        <Route path="/error" element={<Error />} />
        <Route
          exact
          path="/genome-mutation-visualization"
          element={
            <div className="bg-[#f6f7f9] relative">
              <h1 className="text-center pt-4 pb-0 font-bold text-xl">
                {selectedAlgorithm === "PRIEST"
                  ? "PRIEST Spike Site Annotation"
                  : "Genome Sequence Mutation Visualization"}
              </h1>
              <div className="absolute top-0 flex justify-center items-center ">
                <img
                  src={logo}
                  className="w-[7rem] h-auto ml-[5.5rem]"
                  alt="Covidmutext Logo"
                />
              </div>
              {selectedAlgorithm === "PRIEST" ? (
                <SpikeMutationPanel annotation={priestAnnotation} />
              ) : (
                <div className="block md:flex md:justify-normal">
                  {genomeData && genomeData.length > 0 && (
                    <GenomeChart
                      genomeData={genomeData}
                      genomeSequence={genomeSequence}
                    />
                  )}
                  {!selectedProteinRegion && (
                    <DoughnutChart data={protein_mutation_probs} />
                  )}
                </div>
              )}
            </div>
          }
        />
        <Route
          exact
          path="/case-studies/known-hotspot"
          element={<KnownHotspotCaseStudyPage />}
        />
        <Route exact path="/contact-us" element={<Contact />} />
        <Route exact path="/about" element={<About />} />
      </Routes>
    </div>
  );
}

export default App;
