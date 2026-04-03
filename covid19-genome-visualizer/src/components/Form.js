import React from "react";
import { nodeIds } from "../data/nodeIds";

const Form = () => {
  return (
    <div className="p-8 flex justify-center w-full">
      <form className="h-screen">
        <p className="required p-4">{`Prediction Model ${"*"}`}</p>
        <select className="w-full p-4" placeholder="">
          <option disabled selected hidden>
            select prediction model{" "}
          </option>
        </select>

        <select className="w-full p-4">
          <option disabled selected hidden>
            Covid19 Variant Id{" "}
          </option>
          {nodeIds.map((item) => {
            return (
              <>
                <option>{item}</option>
              </>
            );
          })}
        </select>
        <div className="flex justify-between gap-x-5 ">
          <div>
            <label className="p-4">Elasped Days</label>
            <div>
              <input
                type="text"
                className="p-4 border-[2px] border-black"
              ></input>
            </div>
          </div>
          <select className="">
            <option disabled selected>
              Protein Region
            </option>
          </select>
        </div>
      </form>
    </div>
  );
};

export default Form;
