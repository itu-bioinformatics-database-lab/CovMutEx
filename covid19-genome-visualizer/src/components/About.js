import React from "react";
import image1 from "../CovMutex.About-images-0.jpg";
import image2 from "../CovMutex.About-images-1.jpg";
import image3 from "../CovMutex.About-images-2.jpg";
import image4 from "../CovMutex.About-images-3.jpg";
import demo from "../CovMutExVideo.mp4";
import TutorialPage from "./TutorialPage"; 

export const About = () => {
  return (
    // The main container for the page
    <div className="flex justify-center items-center py-8">
      <div className="flex flex-col items-center w-full max-w-5xl px-4">
      
        <h1 className="font-extrabold text-4xl mb-8 text-center"> 
          See CovMutEx in Action 
        </h1>
        
        
        <video src={demo} controls className="w-full rounded-lg shadow-lg" />

        <div className="my-16 w-full flex justify-center">
          <TutorialPage />
        </div>

        
        <h2 className="font-extrabold text-3xl mb-8 text-center">
          Prediction Model Architectures
        </h2>
        <div className="w-full space-y-4">
          <img src={image1} alt="CovMutEx Interface Screenshot 1" className="w-full rounded-lg shadow-lg" />
          <img src={image2} alt="CovMutEx Interface Screenshot 2" className="w-full rounded-lg shadow-lg" />
          <img src={image3} alt="CovMutEx Interface Screenshot 3" className="w-full rounded-lg shadow-lg" />
          <img src={image4} alt="CovMutEx Interface Screenshot 4" className="w-full rounded-lg shadow-lg" />
        </div>

      </div>
    </div>
  );
};
