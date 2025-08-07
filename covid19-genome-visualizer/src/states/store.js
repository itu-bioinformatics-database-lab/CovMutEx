import { configureStore } from "@reduxjs/toolkit";
import genomeReducer from "../features/genome/genomeSlice";

export const store = configureStore({
	reducer: {
		genome: genomeReducer,
	},
});
