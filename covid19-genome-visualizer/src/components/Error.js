import React from "react";


const Error = () => {
  return (
    <div className="min-h-screen bg-blue-50 flex items-center justify-center px-4 py-8">
      <div className="max-w-md w-full space-y-6 text-center">
        <div className="bg-blue-100 rounded-full w-24 h-24 mx-auto flex items-center justify-center">
          <span className="text-6xl font-bold text-blue-600">404</span>
        </div>
        <div>
          <h1 className="text-3xl font-bold text-blue-800 mb-4">
            Page Not Found
          </h1>
          <p className="text-blue-600 mb-6">
            Oops! An error was encountered while processing your request. Please
            try again later.
          </p>
          <p className="text-red-600 mb-6 text-sm font-semibold">
            Error: Unable to retrieve genome data
          </p>
          <a
            href="/"
            className="inline-flex items-center bg-blue-500 text-white px-6 py-3 rounded-lg hover:bg-blue-600 transition duration-300"
          >
            Return Home
          </a>
        </div>
      </div>
    </div>
  );
};

export default Error;
