const fs = require('fs');
const packageJson = JSON.parse(fs.readFileSync('./package.json', 'utf8'));

// Set correct versions for key dependencies
const updatedDependencies = {
  "react": "18.2.0",
  "react-dom": "18.2.0",
  "@material-tailwind/react": "2.1.8",
  "tailwindcss": "3.3.5",
  "postcss": "8.4.31",
  "autoprefixer": "10.4.16",
  "@tailwindcss/postcss": "0.1.1"
};

// Update dependencies
if (packageJson.dependencies) {
  Object.keys(updatedDependencies).forEach(dep => {
    if (packageJson.dependencies[dep]) {
      packageJson.dependencies[dep] = updatedDependencies[dep];
    }
  });
}

// Update devDependencies
if (packageJson.devDependencies) {
  Object.keys(updatedDependencies).forEach(dep => {
    if (packageJson.devDependencies[dep]) {
      packageJson.devDependencies[dep] = updatedDependencies[dep];
    }
  });
}

// Ensure tailwind related packages are in devDependencies
if (!packageJson.devDependencies) {
  packageJson.devDependencies = {};
}

// Move Tailwind dependencies to devDependencies if they're in dependencies
["tailwindcss", "postcss", "autoprefixer", "@tailwindcss/postcss"].forEach(dep => {
  if (packageJson.dependencies && packageJson.dependencies[dep]) {
    packageJson.devDependencies[dep] = packageJson.dependencies[dep];
    delete packageJson.dependencies[dep];
  }
});

// Add resolutions to ensure consistent React version
packageJson.resolutions = {
  "react": "18.2.0",
  "react-dom": "18.2.0"
};

fs.writeFileSync('./package.json', JSON.stringify(packageJson, null, 2));
