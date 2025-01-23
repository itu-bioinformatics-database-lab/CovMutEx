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
  const API_URL = process.env.REACT_APP_API_URL;

  // Reset `selectedProteinRegion` on navigation
  useEffect(() => {
    if (location.pathname === "/") {
      dispatch(resetProteinRegion());
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

  const onSubmit = async (nodeId, elapsedDay, selectedModel, event) => {
    try {
      const response = await fetch(`${API_URL}/api/predict/`, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
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

      setGenomeSequence(data.genomeSequence);
      setGenomeData(data.genomeData);
      setProteinRegionPossibilities(data.proteinRegionPossibilities);
      setProteinMutationProbs(data.protein_mutation_probs);

      dispatch(
        setDataset({
          dataset: data.genomeData,
          genome: data.genomeSequence,
          pr_poss: data.proteinRegionPossibilities,
          isSelected: selectedProteinRegion !== null,
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
    <div className="overflow-y-hidden">
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
                Genome Sequence Mutation Visualization
              </h1>
              <div className="absolute top-0 flex justify-center items-center ">
                <img
                  src={logo}
                  className="w-[7rem] h-auto ml-[5.5rem]"
                  alt="Covidmutext Logo"
                />
              </div>
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
            </div>
          }
        />
        <Route exact path="/contact-us" element={<Contact />} />
        <Route exact path="/about" element={<About />} />
      </Routes>
    </div>
  );
}

export default App;
