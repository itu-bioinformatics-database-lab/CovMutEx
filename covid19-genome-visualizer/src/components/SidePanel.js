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
  onProteinClick, 
}) => {
  return (
    <div className="lg:mr-4 mt-6 w-full lg:w-auto">
      <Card className="w-full lg:w-60 h-auto lg:h-[90vh] shadow-md border-[2px] py-5 overflow-y-auto">
        <div className="flex justify-center items-center">
          <Typography variant="h5" color="blue-gray" className="mb-2">
            Protein Regions
          </Typography>
          <TooltipCustom />
        </div>
        <div className="m-2 w-full flex justify-between px-4">
          <Switch
            onChange={handleShowFullAnnotation}
            label="Show Protein Regions"
            labelProps={{ className: "font-medium !min-w-[40px]" }}
            containerProps={{
              className: "!bg-gray-300",
            }}
            color="blue"
          />
        </div>

        <List className="px-2">
          {Object.keys(proteinRegions).map((protein) => (
            <ListItem
              key={protein}
              // Add cursor style, hover effect, and the onClick handler
              className="text-sm p-1 mb-[0.4rem] flex items-center justify-between font-normal border-solid border-2 cursor-pointer hover:bg-gray-100"
              onClick={() => onProteinClick(protein)}
              onMouseEnter={() => {
                onProteinHover(protein);
              }}
              onMouseLeave={() => {
                onProteinLeave();
              }}
            >
              <span className="font-semibold truncate max-w-[50%]">
                {protein}:
              </span>
              <ListItemSuffix className="ml-4">
                <Chip
                  value={proteinRegions[protein]}
                  variant="ghost"
                  size="sm"
                  className="rounded-full px-2 py-1 truncate max-w-[100px]"
                  style={{
                    backgroundColor: proteinRegionColorMap[protein],
                    overflow: "hidden",
                    textOverflow: "ellipsis",
                    whiteSpace: "nowrap",
                  }}
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
