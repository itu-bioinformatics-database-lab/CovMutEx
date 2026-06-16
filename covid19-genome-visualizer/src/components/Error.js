import React from "react";
import { Link } from "react-router-dom";

const Error = () => {
  return (
    <div className="min-h-[calc(100vh-3.5rem)] bg-gradient-to-br from-gray-50 via-blue-50/30 to-gray-50 dark:from-gray-950 dark:via-gray-900 dark:to-gray-950 flex items-center justify-center px-4 py-8 transition-colors">
      <div className="max-w-md w-full text-center animate-fade-in">
        {/* Error Icon */}
        <div className="bg-red-100 dark:bg-red-900/30 rounded-full w-24 h-24 mx-auto flex items-center justify-center mb-6">
          <span className="text-5xl font-bold text-red-500 dark:text-red-400">!</span>
        </div>

        <h1 className="text-3xl font-bold text-gray-800 dark:text-gray-100 mb-3">
          Something Went Wrong
        </h1>
        <p className="text-gray-500 dark:text-gray-400 mb-2">
          An error was encountered while processing your request.
        </p>
        <p className="text-red-500 dark:text-red-400 mb-8 text-sm font-medium bg-red-50 dark:bg-red-900/20 rounded-xl px-4 py-2 inline-block">
          Unable to retrieve genome data
        </p>

        <div className="flex flex-col sm:flex-row gap-3 justify-center">
          <Link
            to="/"
            className="inline-flex items-center justify-center bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 text-white px-6 py-3 rounded-xl font-semibold shadow-lg shadow-blue-500/20 hover:shadow-blue-500/40 transition-all"
          >
            Return Home
          </Link>
          <button
            onClick={() => window.location.reload()}
            className="inline-flex items-center justify-center bg-gray-200 hover:bg-gray-300 dark:bg-gray-800 dark:hover:bg-gray-700 text-gray-700 dark:text-gray-300 px-6 py-3 rounded-xl font-semibold transition-all"
          >
            Try Again
          </button>
        </div>
      </div>
    </div>
  );
};

export default Error;
