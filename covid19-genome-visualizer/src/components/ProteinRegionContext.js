import React, { createContext, useState } from "react";

export const ProteinRegionContext = createContext();

export const ProteinRegionProvider = ({ children }) => {
  const [selectedProteinRegion, setSelectedProteinRegion] = useState(null);

  return (
    <ProteinRegionContext.Provider
      value={{ selectedProteinRegion, setSelectedProteinRegion }}
    >
      {children}
    </ProteinRegionContext.Provider>
  );
};
