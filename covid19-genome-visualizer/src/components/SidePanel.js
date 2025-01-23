import React from "react";
import {
  List,
  ListItem,
  Card,
  Typography,
  ListItemSuffix,
  Chip,
  Switch,
} from "@material-tailwind/react";
import { proteinRegionColorMap } from "../utils/proteinRegionColorMap";
import { TooltipCustom } from "./TooltipComponent";

const SidePanel = ({
  proteinRegions,
  onProteinHover,
  onProteinLeave,
  handleShowFullAnnotation,
}) => {
  return (
    <div className="mr-4 mt-6">
      <Card className="w-60 h-[90vh] shadow-md border-[2px] py-5">
        <div className="flex justify-center">
          <Typography variant="h5" color="blue-gray" className="mb-2">
            Protein Regions
          </Typography>
          <TooltipCustom />
        </div>
        <div className="m-2 flex justify-center">
          <Switch
            onChange={handleShowFullAnnotation}
            label="Show Protein Regions"
            labelProps={{ className: "font-medium" }}
            color="blue"
          />
        </div>

        <List>
          {Object.keys(proteinRegions).map((protein) => (
            <ListItem
              key={protein}
              className="text-sm p-1.5 mb-[0.4rem] flex items-center justify-end font-normal border-solid border-2"
              onMouseEnter={() => {
                // console.log("Hovering over:", protein); // Debug log
                onProteinHover(protein);
              }}
              onMouseLeave={() => {
                // console.log("Leaving:", protein); // Debug log
                onProteinLeave();
              }}
            >
              <span className="font-semibold">{protein}:</span>
              <ListItemSuffix>
                <Chip
                  value={proteinRegions[protein]}
                  variant="ghost"
                  size="sm"
                  className="rounded-full"
                  style={{ backgroundColor: proteinRegionColorMap[protein] }}
                />
              </ListItemSuffix>
            </ListItem>
          ))}
        </List>
      </Card>
    </div>
  );
};

export default SidePanel;
