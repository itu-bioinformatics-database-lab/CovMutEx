import React, { useState } from "react";
import { nodeIds as nodes } from "../data/nodeIds";
import { modelList } from "../data/modelList";
import { useDispatch, useSelector } from "react-redux";
import { Select as Select2, Button, Input } from "@material-tailwind/react";
import { MdOutlineCreate } from "react-icons/md";
import Select from "react-select";
import LoadingSpinner from "./Spinner";
import { proteinRegions } from "../data/proteinRegions";
import DropDown from "./DropDown";
import { useNavigate } from "react-router-dom";
import {
  updateProteinRegion,
  resetProteinRegion,
} from "../features/genome/genomeSlice";

import logo from "../CovMutexLogo-removebg-preview.png";

const customStyles = {
  placeholder: (provided) => ({
    ...provided,
    color: "#BDC3D4",
    fontSize: "14px",
  }),
  control: (provided) => ({
    ...provided,
    minHeight: "40px",
  }),
};

function Navbar({ onNodeSelect, onSubmit }) {
  const dispatch = useDispatch();
  const node = useSelector((state) => state.genome.nodeList);
  const selectedProteinRegion = useSelector(
    (state) => state.genome.selectedProteinRegion
  );
  const [_nodeId, setNodeId] = useState(null);
  const [_elapsedDay, setElapsedDay] = useState("");
  const [selectedModel, setSelectedModel] = useState(null);
  const [loading, setLoading] = useState(false);
  const navigate = useNavigate();

  const handleElapsedDayChange = (e) => {
    const value = Math.max(0, Number(e.target.value));
    setElapsedDay(value === 0 ? "" : value.toString());
  };

  const handleProteinRegionChange = (opt) => {
    const selectedRegion = opt ? opt.value : null;
    if (selectedRegion == null || selectedProteinRegion === undefined) {
      dispatch(resetProteinRegion()); // Reset Redux state
    } else {
      console.log("Dispatching selected protein region:", selectedRegion);
      dispatch(updateProteinRegion(selectedRegion));
    }
  };

  const handleSubmit = async (e) => {
    e.preventDefault();
    setLoading(true);
    console.log("Form Submission Triggered");
    console.log("Selected Node ID:", _nodeId);
    console.log("Elapsed Days:", _elapsedDay);
    console.log("Selected Model:", selectedModel);

    if (_nodeId && _elapsedDay && selectedModel) {
      if (!selectedProteinRegion) {
        dispatch(resetProteinRegion());
      }
      console.log("onNodeSelect called with params:", {
        _nodeId,
        _elapsedDay,
        selectedModel,
      });
      onNodeSelect(_nodeId, _elapsedDay, selectedModel);
      await onSubmit(_nodeId, _elapsedDay, selectedModel);
      navigate(`/genome-mutation-visualization`);
    } else {
      alert("Please fill in all required fields before submitting.");
    }
  };

  if (loading || !nodes) {
    return <LoadingSpinner />;
  }

  return (
    <div className="relative">
      <div className="absolute left-1/2 top-0 transform -translate-x-1/2 flex justify-center items-center ">
        <img src={logo} className="w-[20rem] h-[16rem] mb-[2rem]" alt="Covidmutext Logo" />
      </div>

      <div className="min-h-screen p-4 flex justify-center items-center ">
        <form
          className="w-full navbar max-w-xl space-y-4"
          onSubmit={handleSubmit}
        >
          {/* Prediction Model */}
          <div className="w-full">
            <label className="text-sm mb-1 text-blue-600 font-semibold block">
              Prediction Model{" "}
              <span className="text-sm text-red-300 font-semibold">*</span>
            </label>
            <Select
              required
              onChange={(opt) => setSelectedModel(opt.value)}
              className="w-full text-sm"
              options={modelList.map((model) => ({
                label: model.name,
                value: model.path,
              }))}
              placeholder="Select Prediction Model"
              styles={customStyles}
              name="model"
            />
          </div>

          {/* Covid19 Variant Id */}
          <div className="w-full">
            <label className="text-sm mb-1 text-blue-600 font-semibold block">
              Covid19 Variant Id{" "}
              <span className="text-sm text-red-300 font-semibold">*</span>
            </label>
            <DropDown items={nodes} setNodeId={setNodeId} />
          </div>

          {/* Elapsed Days and Protein Region Row */}
          <div className="flex flex-col sm:flex-row gap-4">
            {/* Elapsed Days */}
            <div className="flex-1">
              <label className="text-sm mb-1 font-semibold block text-blue-600">
                Elapsed Days{" "}
                <span className="text-sm text-red-300 font-semibold">*</span>
              </label>
              <Input
                required
                value={_elapsedDay}
                onChange={handleElapsedDayChange}
                className="w-full h-10 border-black"
                type="number"
                min={1}
                label="e.g., 120"
              />
            </div>

            {/* Protein Region */}
            <div className="flex-1">
              <label className="text-sm mb-1 text-blue-600 font-semibold block">
                Select Protein Region
              </label>
              <Select
                className="w-full text-sm"
                onChange={handleProteinRegionChange}
                options={Object.keys(proteinRegions).map((pr) => ({
                  label: pr,
                  value: pr,
                }))}
                styles={customStyles}
                placeholder="Optional"
              />
            </div>
          </div>

          {/* Submit Button */}
          <div className="w-full pt-4">
            <Button
              color="blue"
              variant="gradient"
              type="submit"
              className="w-full flex justify-center items-center gap-2"
            >
              Predict <MdOutlineCreate className="h-4 w-4" />
            </Button>
          </div>
        </form>
      </div>
    </div>
  );
}

export default Navbar;

