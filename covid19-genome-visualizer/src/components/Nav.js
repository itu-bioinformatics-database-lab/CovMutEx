import React from "react";
import { Link } from "react-router-dom";

const Nav = () => {
  return (
    <nav className="bg-blue-600 p-4">
      <ul className="flex justify-self-start text-white">
        <li>
          <Link to="/" className="hover:text-blue-200">
            Home
          </Link>
        </li>
        <li>
          <Link to="/benchmark" className="hover:text-blue-200 ml-[3rem]">
            Benchmark
          </Link>
        </li>
        <li>
          <Link to="/about" className="hover:text-blue-200 ml-[3rem]">
            About
          </Link>
        </li>
        <li>
          <Link to="/contact-us" className="hover:text-blue-200 ml-[3rem]">
            Contact
          </Link>
        </li>
      </ul>
    </nav>
  );
};

export default Nav;