import React from "react";
import { Link } from "react-router-dom";

const Nav = () => {
  return (
    <nav className="bg-blue-600 p-4">
      <ul className="flex justify-around text-white">
        <li>
          <Link to="/" className="hover:text-blue-200">
            Home
          </Link>
        </li>
        <li>
          <Link to="/about" className="hover:text-blue-200">
            About
          </Link>
        </li>
        <li>
          <Link to="/contact-us" className="hover:text-blue-200">
            Contact
          </Link>
        </li>
      </ul>
    </nav>
  );
};

export default Nav;
