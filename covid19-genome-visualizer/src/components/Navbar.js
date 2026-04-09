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
  const defaultModel =
    modelList.find((model) => model.path === "balanced_data_model") ||
    modelList[0] ||
    null;
  const [selectedModel, setSelectedModel] = useState(defaultModel?.path || null);
  const [loading, setLoading] = useState(false);
  const navigate = useNavigate();
  const isPriestSelected = selectedModel === "PRIEST";

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

    if (_nodeId && selectedModel && (isPriestSelected || _elapsedDay)) {
      if (isPriestSelected || !selectedProteinRegion) {
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
    <div className="min-h-screen p-4 flex flex-col justify-center items-center">
      {/* Centered Logo */}
      <div className="flex justify-center items-center w-full">
        <img
          src={logo}
          alt="Covidmutext Logo"
          className="w-[18rem] h-[9rem] sm:w-64 sm:h-[9rem] md:w-80 md:h-[10rem] lg:w-[22rem] lg:h-[14rem] xl:w-[32rem] xl:h-[20rem] object-contain"
        />
      </div>

      {/* Form */}
      <form
        className="w-full navbar max-w-xl space-y-4"
        onSubmit={handleSubmit}
      >
        {/* Prediction Model */}
        <div className="w-full">
          <label className="text-sm mb-1 text-blue-600 font-semibold block">
            Algorithm{" "}
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
            
            defaultValue = { defaultModel ? {
            label : defaultModel.name,
            value: defaultModel.path
            
            } : null
            }
            placeholder="Select Algorithm"
            styles={customStyles}
            name="model"
          />
        </div>

        {/* Covid19 Variant Id */}
        <div className="w-full">
          <label className="text-sm mb-1 text-blue-600 font-semibold block focus:outline-none focus:ring-2 focus:ring-blue-400 focus:border-blue-400">
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
              {!isPriestSelected ? (
                <span className="text-sm text-red-300 font-semibold">*</span>
              ) : null}
            </label>
            <Input
              required={!isPriestSelected}
              value={_elapsedDay}
              onChange={handleElapsedDayChange}
              type="number"
              min={1}
              label="e.g., 120"
              className="!border !border-gray-300 !rounded-md focus:!border-blue-500 focus:!border-2 !text-sm !font-normal focus:!ring-0"
              containerProps={{
                className: "min-h-[40px]",
              }}
            />
          </div>

          {/* Protein Region */}
          <div className="flex-1">
            <label className="text-sm mb-1 text-blue-600 font-semibold block">
              Select Protein Region
            </label>
            <Select
              className="w-full text-sm"
              isDisabled={isPriestSelected}
              onChange={handleProteinRegionChange}
              options={[ { label: "\u00A0", value: "" }, ...Object.keys(proteinRegions).map((pr) => ({ label: pr, value: pr, })) ]}
              styles={customStyles}
              placeholder={isPriestSelected ? "Not used for PRIEST" : "Optional"}
            />
          </div>
        </div>

        {/* Submit Button */}
        <div className="w-full pt-4">
          <Button
            color="blue"
            variant="gradient"
            type="submit"
            className="w-full bg-blue-500 flex justify-center items-center gap-2"
          >
            Explore Mutations <MdOutlineCreate className="h-4 w-4" />
          </Button>
        </div>
      </form>
    </div>
  );
}

export default Navbar;
