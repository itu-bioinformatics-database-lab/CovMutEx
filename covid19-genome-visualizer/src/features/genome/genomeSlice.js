import { createSlice, current } from "@reduxjs/toolkit";
import { generateRandomSequence } from "../../helpers/fooGenomeElements";
import { proteinRegions } from "../../data/proteinRegions";
import axios from "axios";

const initialState = {
  randomSeq: null,
  possibilityMap: null,
  windowSlices: null,
  proteinRegionPossMap: null,
  chartData: null,
  chartTitle: "Full Sequence",
  isWholeSequenceSelected: true,
  elapsedDay: 0,
  nodeId: "",
  model: "",
  modelList: [],
  nodeList: [],
  seq: null,
  realChartData: null,
  proteinRegionPossMap2: null,
  showDoughnut: false,
  selectedProteinRegion: null,
  loading: false,
  dataset: [],
  mutationProbabilities: [],
  genome: "",
  pr_poss: {},
  isSelected: false,
};

export const genomeSlice = createSlice({
  name: "genome",
  initialState,
  reducers: {
    generate: (state) => {
      state.randomSeq = generateRandomSequence();
    },
    setDataset: (state, action) => {
      
      const {
        dataset = [],
        genome = "",
        pr_poss = {},
        isSelected = false,
        selectedProteinRegion = null,
      } = action.payload;

      
      state.dataset = [...dataset];
      state.genome = genome;
      state.pr_poss = { ...pr_poss };
      state.isSelected = isSelected;
      state.selectedProteinRegion = selectedProteinRegion;
      state.showDoughnut = selectedProteinRegion === null;
      state.chartData = [...dataset];
      state.realChartData = [...dataset];
      
    },

    showProteinRegion: (state, action) => {
      
      const region = proteinRegions[action.payload];
      if (!region) {
        console.error(`Protein region "${action.payload}" not found.`);
        return;
      }
      const [start, end] = region.split("-").map((pos) => parseInt(pos));
      const chartData = current(state).chartData || [];
      state.chartTitle = action.payload;
      state.realChartData = chartData.map((arr) => arr.slice(start, end + 1));
      state.isWholeSequenceSelected = false;
      state.selectedProteinRegion = action.payload;
      
    },
    selectNode: (state, action) => {
      
      const [node, elapsedDay, model] = action.payload;
      state.nodeId = node;
      state.elapsedDay = elapsedDay;
      state.model = model;
    },
    resetChart: (state) => {
      
      state.realChartData = [...(current(state).chartData || [])];
      state.isWholeSequenceSelected = true;
      state.chartTitle = "Full Sequence";
      state.selectedProteinRegion = null;
      state.showDoughnut = true;
     
    },
    loadNodesAndModels: (state, action) => {
     
      const [models, nodes] = action.payload;
      state.modelList = [...(models || [])];
      state.nodeList = [...(nodes || [])];
    },
    setLoading: (state, action) => {
      
      state.loading = action.payload;
    },
    updateProteinRegion: (state, action) => {
      
      state.selectedProteinRegion = action.payload || null; // Reset to null if payload is null/undefined
      state.showDoughnut = action.payload === null;
      state.isWholeSequenceSelected = action.payload === null;
      state.chartTitle =
        action.payload === null ? "Full Sequence" : action.payload;
     
    },
    resetProteinRegion: (state) => {
      
      state.selectedProteinRegion = null;
      
    },
  },
});

export const {
  generate,
  setDataset,
  showProteinRegion,
  resetChart,
  selectNode,
  loadNodesAndModels,
  setLoading,
  updateProteinRegion,
  resetProteinRegion,
} = genomeSlice.actions;

export const submitForm =
  (nodeId, elapsedDay, selectedModel, selectedProteinRegion) =>
  async (dispatch) => {
    try {
     
      dispatch(setLoading(true));
      const API_URL = process.env.REACT_APP_API_URL;

      const response = await axios.post(`${API_URL}/api/predict/`, {
        nodeId,
        elapsedDay,
        selectedModel,
        selectedProteinRegion,
      });

      const { dataset, genome, pr_poss } = response.data;
      const isSelected = Boolean(selectedProteinRegion);
      

      

      dispatch(
        setDataset({
          dataset,
          genome,
          pr_poss,
          isSelected,
          selectedProteinRegion,
        })
      );
    } catch (error) {
      console.error("Error in form submission:", error);
    } finally {
      dispatch(setLoading(false));
    }
  };

export default genomeSlice.reducer;

