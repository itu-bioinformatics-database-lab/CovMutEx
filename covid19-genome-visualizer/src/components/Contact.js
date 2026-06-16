import React, { useState } from "react";
import logo from "../CovMutexLogo-removebg-preview.png";
import { MdSend, MdEmail, MdPerson, MdSubject } from "react-icons/md";

const Contact = () => {
  const [formData, setFormData] = useState({
    name: "",
    email: "",
    subject: "",
    message: "",
  });
  const [submitted, setSubmitted] = useState(false);

  const handleChange = (e) => {
    const { name, value } = e.target;
    setFormData({ ...formData, [name]: value });
  };

  const handleSubmit = (e) => {
    e.preventDefault();
    setSubmitted(true);
    setTimeout(() => setSubmitted(false), 3000);
  };

  return (
    <div className="min-h-[calc(100vh-3.5rem)] bg-gradient-to-br from-gray-50 via-blue-50/30 to-gray-50 dark:from-gray-950 dark:via-gray-900 dark:to-gray-950 flex items-center justify-center p-4 transition-colors">
      <div className="w-full max-w-lg animate-fade-in">
        {/* Logo */}
        <div className="flex justify-center mb-6">
          <img src={logo} alt="CovMutEx" className="h-20 w-auto object-contain" />
        </div>

        {/* Card */}
        <div className="bg-white dark:bg-gray-900 rounded-2xl shadow-xl dark:shadow-gray-900/50 border border-gray-200 dark:border-gray-800 p-8">
          <h2 className="text-2xl font-bold text-center text-gray-800 dark:text-gray-100 mb-1">
            Get in Touch
          </h2>
          <p className="text-sm text-center text-gray-500 dark:text-gray-400 mb-6">
            Have questions about CovMutEx? We'd love to hear from you.
          </p>

          {/* Success Message */}
          {submitted && (
            <div className="mb-4 p-3 bg-green-50 dark:bg-green-900/30 border border-green-200 dark:border-green-800 rounded-xl text-green-700 dark:text-green-400 text-sm text-center font-medium animate-slide-up">
              Message sent successfully!
            </div>
          )}

          <form
            action="https://formspree.io/f/xlddapzb"
            method="POST"
            onSubmit={handleSubmit}
            className="space-y-4"
          >
            {/* Name */}
            <div>
              <label htmlFor="name" className="text-xs font-bold text-gray-500 dark:text-gray-400 uppercase mb-1 block">
                Your Name
              </label>
              <div className="relative">
                <MdPerson className="absolute left-3 top-1/2 -translate-y-1/2 text-gray-400" />
                <input
                  type="text"
                  id="name"
                  name="name"
                  value={formData.name}
                  onChange={handleChange}
                  className="w-full pl-10 pr-4 py-2.5 border-2 border-gray-200 dark:border-gray-700 rounded-xl bg-white dark:bg-gray-800 text-gray-800 dark:text-gray-200 focus:border-blue-500 focus:ring-2 focus:ring-blue-500/20 transition-all text-sm"
                  placeholder="John Doe"
                  required
                />
              </div>
            </div>

            {/* Email */}
            <div>
              <label htmlFor="email" className="text-xs font-bold text-gray-500 dark:text-gray-400 uppercase mb-1 block">
                Email Address
              </label>
              <div className="relative">
                <MdEmail className="absolute left-3 top-1/2 -translate-y-1/2 text-gray-400" />
                <input
                  type="email"
                  id="email"
                  name="email"
                  value={formData.email}
                  onChange={handleChange}
                  className="w-full pl-10 pr-4 py-2.5 border-2 border-gray-200 dark:border-gray-700 rounded-xl bg-white dark:bg-gray-800 text-gray-800 dark:text-gray-200 focus:border-blue-500 focus:ring-2 focus:ring-blue-500/20 transition-all text-sm"
                  placeholder="you@example.com"
                  required
                />
              </div>
            </div>

            {/* Subject */}
            <div>
              <label htmlFor="subject" className="text-xs font-bold text-gray-500 dark:text-gray-400 uppercase mb-1 block">
                Subject
              </label>
              <div className="relative">
                <MdSubject className="absolute left-3 top-1/2 -translate-y-1/2 text-gray-400" />
                <input
                  type="text"
                  id="subject"
                  name="subject"
                  value={formData.subject}
                  onChange={handleChange}
                  className="w-full pl-10 pr-4 py-2.5 border-2 border-gray-200 dark:border-gray-700 rounded-xl bg-white dark:bg-gray-800 text-gray-800 dark:text-gray-200 focus:border-blue-500 focus:ring-2 focus:ring-blue-500/20 transition-all text-sm"
                  placeholder="What's this about?"
                  required
                />
              </div>
            </div>

            {/* Message */}
            <div>
              <label htmlFor="message" className="text-xs font-bold text-gray-500 dark:text-gray-400 uppercase mb-1 block">
                Message
              </label>
              <textarea
                id="message"
                name="message"
                value={formData.message}
                onChange={handleChange}
                className="w-full px-4 py-2.5 border-2 border-gray-200 dark:border-gray-700 rounded-xl bg-white dark:bg-gray-800 text-gray-800 dark:text-gray-200 focus:border-blue-500 focus:ring-2 focus:ring-blue-500/20 transition-all text-sm resize-none"
                rows="4"
                placeholder="Tell us more..."
                required
              />
            </div>

            {/* Submit */}
            <button
              type="submit"
              className="w-full bg-blue-600 hover:bg-blue-700 dark:bg-blue-500 dark:hover:bg-blue-600 text-white py-3 rounded-xl font-semibold flex items-center justify-center gap-2 shadow-lg shadow-blue-500/20 hover:shadow-blue-500/40 transition-all"
            >
              <MdSend className="text-lg" />
              Send Message
            </button>
          </form>
        </div>

        {/* Footer hint */}
        <p className="text-center text-xs text-gray-400 dark:text-gray-600 mt-4">
          ITU Bioinformatics & Database Lab
        </p>
      </div>
    </div>
  );
};

export default Contact;
