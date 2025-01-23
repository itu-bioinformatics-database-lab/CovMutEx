import React from "react";
import image1 from "../CovMutex.About-images-0.jpg";
import image2 from "../CovMutex.About-images-1.jpg";
import image3 from "../CovMutex.About-images-2.jpg";
import image4 from "../CovMutex.About-images-3.jpg";

export const About = () => {
  return (
    <div className="flex justify-center items-center">
      <div className="flex flex-col items-center">
        <img src={image1} alt="about" className="w-[70rem]" />
        <img src={image2} alt="about" className="w-[70rem]" />
        <img src={image3} alt="about" className="w-[70rem]" />
        <img src={image4} alt="about" className="w-[70rem]" />
      </div>
    </div>
  );
};
