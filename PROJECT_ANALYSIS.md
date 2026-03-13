# CovMutEx - Comprehensive Project Analysis

## Executive Summary

**CovMutEx** (COVID-19 Mutation Explorer) is a full-stack web application designed to predict and visualize future SARS-CoV-2 mutations. The tool combines machine learning-based mutation prediction with interactive genome visualization, enabling researchers to identify mutation "hotspots" and analyze protein regions that are more prone to mutations.

---

## 🏗️ Project Architecture Overview

### System Architecture Type: **Client-Server with REST API**

```
┌─────────────────────────────────────────────────────────────┐
│                      USER INTERFACE                          │
│              (React SPA - Port 3000)                         │
│  - Interactive Charts (Chart.js, Recharts)                   │
│  - Protein Region Visualization                              │
│  - Mutation Probability Display                              │
└────────────────┬────────────────────────────────────────────┘
                 │ HTTP/REST API
                 │ (JSON)
┌────────────────▼────────────────────────────────────────────┐
│                   BACKEND API SERVER                         │
│              (Django + DRF - Port 8000)                      │
│  - Feature Extraction Engine                                 │
│  - TensorFlow Model Inference                                │
│  - Phylogenetic Analysis                                     │
│  - Genome Data Processing                                    │
└────────────────┬────────────────────────────────────────────┘
                 │
        ┌────────┴────────┐
        │                 │
┌───────▼──────┐  ┌──────▼──────────┐
│  ML Models   │  │  Data Files     │
│  (.keras)    │  │  - genome.txt   │
│              │  │  - mutations.txt│
│              │  │  - phylo tree   │
└──────────────┘  └─────────────────┘
```

---

## 📦 Project Structure

### Root Level
```
CovMutEx/
├── README.md                          # Project overview and description
├── covid19-genome-visualizer/         # Frontend React application
└── genome_extractor/                  # Backend Django application
```

---

## 🎨 Frontend Application (`covid19-genome-visualizer/`)

### Technology Stack
- **Framework:** React 18.2.0
- **State Management:** Redux Toolkit (@reduxjs/toolkit)
- **UI Library:** Material Tailwind React
- **Styling:** Tailwind CSS 4.1.7 with PostCSS
- **Routing:** React Router DOM v6.14.2
- **Charting Libraries:**
  - Chart.js 4.4.0 (primary)
  - Recharts 2.10.3
  - Plotly.js 2.27.1
  - ECharts 5.4.3
  - D3.js 7.8.5
- **HTTP Client:** Axios 1.6.3
- **Build Tool:** React Scripts 5.0.1 (Create React App)

### Architecture Pattern
**Flux Architecture with Redux**
- Unidirectional data flow
- Centralized state management
- Action-based state mutations

---

### File-by-File Breakdown

#### Entry Points
- **`public/index.html`**: HTML template with root div
- **`src/index.js`**: 
  - React application entry point
  - Wraps App with Redux Provider, React Router, and Material Tailwind ThemeProvider
  - Renders to DOM root element

#### Core Application Files

**`src/App.js`** (Main Application Component)
- **Purpose:** Root component managing routing and application state
- **Key Responsibilities:**
  - Route configuration (home, visualization, error, contact, about)
  - API communication with backend (`/api/predict/`)
  - Handles dual-fetch strategy:
    1. First fetch: Gets protein-specific calculations
    2. Second fetch (conditional): Gets full genome data for chart rendering
  - Manages loading states and navigation
  - Prevents accidental page refreshes when data is loaded
- **Key Features:**
  - Protein region selection and reset
  - Form submission handling (`onSubmit` function)
  - Data normalization before passing to visualization
  - Error boundary handling with navigation to error page

**`src/config.js`**
- **Purpose:** Configuration file for API endpoints
- **Content:** API base URL (`http://127.0.0.1:8000/`)
- **Usage:** Can be overridden by environment variable `REACT_APP_API_URL`

#### State Management (`src/states/` and `src/features/`)

**`src/states/store.js`**
- **Purpose:** Redux store configuration
- **Structure:** 
  - Single reducer: `genome` (from genomeSlice)
  - Configured with Redux Toolkit's `configureStore`

**`src/features/genome/genomeSlice.js`** (Core State Management)
- **Purpose:** Redux slice managing all genome-related state
- **State Structure:**
  ```javascript
  {
    randomSeq: null,              // Random sequence generator
    possibilityMap: null,          // Mutation possibilities
    windowSlices: null,            // Genome window data
    chartData: null,               // Full chart dataset
    realChartData: null,           // Filtered/processed chart data
    chartTitle: string,            // Display title
    isWholeSequenceSelected: bool, // Full vs. protein region
    elapsedDay: number,            // Days elapsed for prediction
    nodeId: string,                // Selected phylogenetic node
    model: string,                 // Selected ML model path
    modelList: array,              // Available models
    nodeList: array,               // Available phylogenetic nodes
    seq: null,                     // Genome sequence
    selectedProteinRegion: string, // Currently selected protein
    showDoughnut: bool,            // Show/hide doughnut chart
    loading: bool,                 // Loading state
    dataset: array,                // Raw prediction dataset
    genome: string,                // Full genome sequence
    pr_poss: object,               // Protein region possibilities
    isSelected: bool               // Whether user has made selections
  }
  ```

- **Key Actions/Reducers:**
  - `generate()`: Generate random sequence
  - `setDataset()`: Store prediction results from backend
  - `showProteinRegion(proteinName)`: Filter data to specific protein region
  - `resetChart()`: Reset to full genome view
  - `selectNode([nodeId, elapsedDay, model])`: Store user selections
  - `loadNodesAndModels([models, nodes])`: Populate dropdown options
  - `updateProteinRegion(region)`: Update selected protein region
  - `resetProteinRegion()`: Clear protein region selection
  
- **Async Thunk:**
  - `submitForm()`: Handles API call to backend (legacy, mostly handled in App.js now)

**`src/features/genome/genome.js`**
- **Purpose:** Additional genome-related logic and selectors
- **Content:** Helper functions for genome data manipulation

---

#### Component Architecture

**Navigation Components**

1. **`src/components/Nav.js`**
   - Top-level navigation bar
   - Links to Home, Visualization, About, Contact pages
   - Logo display

2. **`src/components/Navbar.js`** (Main Form Component)
   - **Purpose:** User input form for prediction parameters
   - **Key Features:**
     - Node ID selection (dropdown with search/virtualization)
     - Elapsed days input (numeric, minimum 0)
     - Model selection dropdown
     - Protein region filter (optional)
     - Form validation
     - Loading spinner during submission
   - **Technologies:**
     - Material Tailwind components
     - React Select (virtualized for performance with large lists)
     - Redux integration for state updates
   - **Data Flow:**
     ```
     User Input → Form Validation → Redux State Update → 
     API Call (via App.js) → Navigate to Visualization
     ```

**Visualization Components**

3. **`src/components/Recharts.js`** (Main Genome Chart Component)
   - **Purpose:** Primary interactive genome visualization
   - **Chart Type:** Bar chart with extensive customization
   - **Key Features:**
     - Dynamic zooming and panning (chartjs-plugin-zoom)
     - Protein region annotations (colored background zones)
     - Nucleotide-specific coloring
     - Position-level detail on zoom
     - High-resolution view below 1000 nucleotides
     - WebLogo generation for detailed views
     - Custom plugins:
       - `highResLoadingOverlay`: Shows loading state
       - `weblogoOverlay`: Displays WebLogo images
       - `chartVisibilityController`: Manages chart/logo visibility
   - **State Management:**
     - Local state for active protein, zoom levels, view modes
     - Redux selectors for genome data
   - **Interactivity:**
     - Click on bars to zoom
     - Hover for tooltips
     - Protein region highlighting on hover
     - Switch between bar chart and WebLogo views
   - **Performance Optimizations:**
     - Conditional rendering based on zoom level
     - Debounced updates
     - Canvas-based rendering (Chart.js)
   - **Data Structure:**
     ```javascript
     genomeData: [
       [P(A), P(T), P(G), P(C)],  // Position 0
       [P(A), P(T), P(G), P(C)],  // Position 1
       // ... for all 29,904 positions
     ]
     ```

4. **`src/components/DoughnutChart.js`**
   - **Purpose:** Protein region mutation probability summary
   - **Chart Type:** Doughnut chart with interactive legend
   - **Key Features:**
     - Shows percentage of mutations per protein region
     - Normalized view (by protein region length)
     - Click to filter genome chart to specific protein
     - Color-coded by protein (consistent with genome chart)
     - Data table with percentages
   - **Toggle Options:**
     - Raw mutation probability
     - Normalized by protein length
   - **Technologies:**
     - React-Chartjs-2 wrapper
     - Material Tailwind Switch component

5. **`src/components/Test.js` / `ZoomChart`**
   - **Purpose:** Alternative zoom-focused visualization component
   - **Note:** Likely experimental or backup implementation

6. **`src/components/BarChart.js` / `BarChart2.js`**
   - **Purpose:** Additional bar chart implementations
   - **Note:** May be legacy or alternative visualization approaches

7. **`src/components/SidePanel.js`**
   - **Purpose:** Protein region selector panel
   - **Key Features:**
     - List of all protein regions with color coding
     - Click to zoom to specific protein
     - Hover to highlight protein on main chart
     - Toggle to show/hide protein region annotations
     - Displays protein region ranges (e.g., "266-21555")
   - **Interactivity:**
     - `onProteinClick`: Zooms chart to protein region
     - `onProteinHover`: Highlights protein on chart
     - `onProteinLeave`: Removes highlight
     - `handleShowFullAnnotation`: Toggle annotation visibility

8. **`src/components/TooltipComponent.js`**
   - **Purpose:** Custom tooltip displays
   - **Usage:** Information tooltips for UI elements

**UI Helper Components**

9. **`src/components/ProgressBar.js`**
   - **Purpose:** Loading progress indication
   - **Usage:** During data fetching or processing

10. **`src/components/Spinner.js`**
    - **Purpose:** Loading spinner animation
    - **Usage:** Full-screen loading states

11. **`src/components/ZoomSlider.js`**
    - **Purpose:** Manual zoom control slider
    - **Usage:** Alternative zoom interface for charts

12. **`src/components/DropDown.js`**
    - **Purpose:** Reusable dropdown component
    - **Usage:** Generic dropdown for various selections

**Form and Input Components**

13. **`src/components/Form.js`**
    - **Purpose:** Generic form component
    - **Usage:** Reusable form structure

14. **`src/components/RandomSequence.js`**
    - **Purpose:** Random genome sequence generator
    - **Usage:** Testing or demonstration purposes

**Information Pages**

15. **`src/components/About.js`**
    - **Purpose:** About page with project information
    - **Content:**
      - Video demonstration (CovMutExVideo.mp4)
      - Model architecture diagrams (4 images)
      - Tutorial section
    - **Images:**
      - `CovMutex.About-images-0.jpg` through `CovMutex.About-images-3.jpg`
      - Show model architectures and prediction approaches

16. **`src/components/TutorialPage.js`**
    - **Purpose:** User tutorial/guide
    - **Content:** Step-by-step usage instructions

17. **`src/components/Contact.js`**
    - **Purpose:** Contact information page
    - **Content:** Contact details for project maintainers

18. **`src/components/Error.js`**
    - **Purpose:** Error page component
    - **Usage:** Displayed when API errors occur or invalid routes

**Context Providers**

19. **`src/components/ProteinRegionContext.js`**
    - **Purpose:** React Context for protein region data
    - **Usage:** Provides protein region information across components
    - **Note:** May be redundant with Redux implementation

---

#### Data Files

**`src/data/`**

20. **`proteinRegions.js`**
    - **Purpose:** Protein region definitions
    - **Content:** 
      ```javascript
      {
        ORF1ab: "266-21555",
        S: "21563-25384",
        ORF3a: "25393-26220",
        E: "26245-26472",
        M: "26523-27191",
        ORF6: "27202-27387",
        ORF7a: "27394-27759",
        ORF7b: "27756-27887",
        ORF8: "27894-28259",
        N: "28274-29533",
        ORF10: "29558-29674"
      }
      ```
    - **Exports:**
      - `proteinRegions`: Position ranges
      - `proteinRegionsSize`: Calculated lengths

21. **`nodeIds.js`**
    - **Purpose:** List of phylogenetic node IDs
    - **Usage:** Populate node selection dropdown
    - **Content:** Array of node ID strings from phylogenetic tree

22. **`modelList.js`**
    - **Purpose:** Available ML models
    - **Structure:** 
      ```javascript
      [
        { name: "Model Name", path: "model_filename" },
        // ...
      ]
      ```

23. **`possibilityMap.js`**
    - **Purpose:** Pre-computed mutation possibility mappings
    - **Usage:** Quick lookup for mutation probabilities

---

#### Helper Functions and Utilities

**`src/helpers/`**

24. **`helperFunctions.js`**
    - **Purpose:** General utility functions
    - **Functions:**
      - `getColorForNucleotide(nucleotide)`: Returns color for A/T/G/C
      - `nucleotides`: Array constant `['A', 'T', 'G', 'C']`

25. **`fooGenomeElements.js`**
    - **Purpose:** Genome element generation utilities
    - **Functions:**
      - `generateRandomSequence()`: Creates random genome sequences

**`src/utils/`**

26. **`proteinRegionColorMap.js`**
    - **Purpose:** Color mapping for protein regions
    - **Exports:**
      - `proteinRegionColorMap`: Object mapping protein names to colors
      - `proteinRegionColorMapAnnotations`: Chart.js annotation format
    - **Usage:** Consistent coloring across all visualizations

**`src/constants/`**

27. **`constantVars.js`**
    - **Purpose:** Application-wide constants
    - **Content:** Configuration values, thresholds, etc.

---

#### Configuration Files

28. **`tailwind.config.js`**
    - **Purpose:** Tailwind CSS configuration
    - **Customization:** Color palette, breakpoints, theme extensions

29. **`postcss.config.js`**
    - **Purpose:** PostCSS configuration
    - **Plugins:** Tailwind CSS processor

30. **`package.json.backup`**
    - **Purpose:** Backup of dependencies
    - **Note:** Actual `package.json` may be in root or gitignored
    - **Key Dependencies:**
      - React 18.2.0
      - Redux Toolkit 1.9.7
      - Chart.js 4.4.0 + plugins
      - Axios 1.6.3
      - Material Tailwind 2.1.10

---

### Frontend Data Flow

```
User Interaction (Navbar Form)
        ↓
Redux Action Dispatch (selectNode)
        ↓
API Call to Backend (/api/predict/)
        ↓
Response Processing (App.js)
        ↓
Redux Store Update (setDataset)
        ↓
React Component Re-render
        ↓
Chart Visualization Update (Recharts.js)
```

---

## 🔧 Backend Application (`genome_extractor/`)

### Technology Stack
- **Framework:** Django 5.1.2
- **API:** Django REST Framework 3.15.2
- **ML Framework:** TensorFlow 2.17.0 / Keras 3.6.0
- **Scientific Computing:**
  - NumPy 1.26.4
  - Pandas 2.2.3
  - SciKit-Learn 1.5.2
  - Biopython 1.84
- **Database:** SQLite (db.sqlite3)
- **Caching:** Django Cache Framework + HDF5
- **Visualization:** Matplotlib, Logomaker
- **CORS:** django-cors-headers 4.5.0

### Architecture Pattern
**MVC (Model-View-Controller) / Django MVT**
- Models: Data structures
- Views: Business logic and API endpoints
- Templates: Not used (API-only backend)

---

### File-by-File Breakdown

#### Django Project Configuration

**`manage.py`**
- **Purpose:** Django command-line utility
- **Usage:** 
  - `python manage.py runserver` - Start development server
  - `python manage.py migrate` - Apply database migrations
  - `python manage.py createsuperuser` - Create admin user

**`genome_extractor/` (Project Package)**

31. **`genome_extractor/settings.py`**
    - **Purpose:** Django configuration
    - **Key Settings:**
      - `DEBUG = True` (development mode)
      - `ALLOWED_HOSTS`: Includes production domain `covidmutex.itu.edu.tr`
      - `INSTALLED_APPS`:
        - `django.contrib.*` (admin, auth, sessions, etc.)
        - `rest_framework`
        - `genome` (custom app)
        - `corsheaders`
      - `MIDDLEWARE`:
        - CORS middleware (must be first)
        - Security, sessions, CSRF, auth
      - `DATABASES`: SQLite configuration
      - `CORS_ALLOWED_ORIGINS`: Allows localhost:3000 and production domain
      - `STATIC_URL` and `STATICFILES_DIRS`
    - **Security Notes:**
      - Secret key exposed (should be in environment variables)
      - Debug mode enabled (should be False in production)

32. **`genome_extractor/urls.py`**
    - **Purpose:** Root URL configuration
    - **Routes:**
      - `/admin/` → Django admin interface
      - `/` → Includes `genome.urls`

33. **`genome_extractor/wsgi.py`**
    - **Purpose:** WSGI application for deployment
    - **Usage:** Production servers (Gunicorn, uWSGI)

34. **`genome_extractor/asgi.py`**
    - **Purpose:** ASGI application for async deployment
    - **Usage:** Async servers (Daphne, Uvicorn)

---

#### Main Application (`genome/`)

**Core Django Files**

35. **`genome/apps.py`**
    - **Purpose:** App configuration
    - **Content:** GenomeConfig class

36. **`genome/models.py`**
    - **Purpose:** Database models
    - **Content:** Currently empty (no database models defined)
    - **Note:** Application operates mostly with file-based data

37. **`genome/admin.py`**
    - **Purpose:** Django admin configuration
    - **Content:** Currently empty

38. **`genome/tests.py`**
    - **Purpose:** Unit tests
    - **Status:** No tests implemented yet

39. **`genome/urls.py`**
    - **Purpose:** URL routing for genome app
    - **Routes:**
      - `POST /api/predict/` → `predict_genome` (main prediction endpoint)
      - `POST /generate-weblogo/` → `generate_weblogo` (WebLogo image generation)
      - `GET /` → `home` (fallback)

---

**Core Business Logic Files**

40. **`genome/views.py`** (542 lines - Main API Logic)
    - **Purpose:** API endpoints and request handling
    - **Key Functions:**

    **1. `predict_genome(request)` - Main API Endpoint**
    - **HTTP Methods:** GET, POST
    - **Endpoint:** `/api/predict/`
    - **Input Parameters (JSON):**
      ```json
      {
        "nodeId": "phylogenetic_node_id",
        "elapsedDay": 30,
        "selectedModel": "model_name",
        "selectedProteinRegion": "S" or null
      }
      ```
    - **Process Flow:**
      1. Load genome sequence from `genome.txt`
      2. Parse mutations for selected node from `mutations.txt`
      3. Construct variant genome by applying mutations
      4. Extract features using `cache_precomputed_features()`
      5. Load appropriate ML model (.keras file)
      6. Run predictions (handles both single and multi-input models)
      7. Calculate probability distributions using Ti/Tv ratio
      8. Generate genome data array (29,904 × 4 matrix)
      9. Calculate protein region probabilities
      10. Return JSON response
    - **Output Structure:**
      ```json
      {
        "nodeId": "...",
        "elapsedDay": 30,
        "selectedModel": "...",
        "selectedProteinRegion": "S",
        "genomeSequence": "ATGC...",
        "genomeData": [[P(A), P(T), P(G), P(C)], ...],
        "protein_mutation_probs": {
          "ORF1ab": 0.023,
          "S": 0.045,
          ...
        },
        "proteinRegionPossibilities": {...},
        "model_metadata": {...}
      }
      ```

    **2. `generate_weblogo(request)` - WebLogo Generation**
    - **HTTP Method:** POST
    - **Purpose:** Generate sequence logo visualizations
    - **Input:**
      ```json
      {
        "start": 1,
        "end": 25,
        "probability_matrix": [[...], ...],
        "reference_sequence": "ATGC...",
        "nucleotide_order": ["A", "T", "G", "C"]
      }
      ```
    - **Process:**
      1. Validate probability matrix (Nx4 shape)
      2. Create DataFrame with nucleotide columns
      3. Generate logo using logomaker library
      4. Add position labels with reference nucleotides
      5. Style adjustments (rotation, colors, title)
      6. Save to PNG and return as HTTP response
    - **Output:** PNG image (image/png)

    **3. `_get_biologically_distributed_probs(p_no_mutation, p_mutation, ref_nuc, ti_tv_ratio=2.0)`**
    - **Purpose:** Calculate biologically realistic mutation probabilities
    - **Biological Model:**
      - Transitions (purine↔purine or pyrimidine↔pyrimidine): More likely
      - Transversions (purine↔pyrimidine): Less likely
      - Ti/Tv ratio: 2.0 (configurable)
    - **Algorithm:**
      1. Assign no-mutation probability to reference nucleotide
      2. Calculate transition probability: `prob_ti = prob_tv * ti_tv_ratio`
      3. Calculate transversion probability: `prob_tv = p_mutation / (ti_tv_ratio + 2)`
      4. Distribute probabilities accordingly
      5. Normalize to ensure sum = 1.0
    - **Returns:** Dict with probabilities for A, T, G, C

    **4. `predict_mutations(cache_path, genome_seq, mutations, ...)`**
    - **Purpose:** Coordinate feature extraction and model prediction
    - **Process:**
      1. Call `cache_precomputed_features()` to get feature matrix
      2. Load TensorFlow model
      3. Check model type (single-input vs multi-input)
      4. Run batch prediction
      5. Return prediction array (29,904 × 1 or 29,904 × 2)
    - **Performance Notes:**
      - Uses precomputed features from HDF5 cache
      - Batch prediction for efficiency
      - Handles both single-output and dual-output models

    **5. `calculate_genome_data(genome_seq, position_predictions, selected_protein_region)`**
    - **Purpose:** Convert raw predictions to probability distributions
    - **Logic:**
      - For each genome position:
        - Extract prediction (1 or 2 values)
        - Convert to no-mutation and mutation probabilities
        - Apply biological distribution (Ti/Tv)
        - Store probabilities for A, T, G, C
    - **Protein Region Handling:**
      - If region selected: Only process positions in that region
      - If no region: Process entire genome
    - **Returns:** Array of shape [4, N] where N = genome length or region length

    **6. `calculate_protein_region_probabilities(position_predictions, protein_regions, genome_seq_length)`**
    - **Purpose:** Aggregate mutation probabilities by protein region
    - **Method:**
      - For each protein region:
        - Extract predictions for positions in that region
        - Calculate average probability
        - Store as float (for JSON serialization)
    - **Returns:** Dict mapping protein names to average probabilities

    **Helper Functions:**
    - `measure_time(label, start_time)`: Performance logging
    - `read_genome_sequence(file_path)`: Read and clean genome file
    - `home(request)`: Fallback handler (redirects to predict_genome)

---

41. **`genome/feature_extractor.py`** (1,078 lines - Feature Engineering Engine)
    - **Purpose:** Core feature extraction and preprocessing
    - **Key Functions:**

    **Feature Extraction Functions:**

    **1. `get_sequence()`**
    - **Purpose:** Load reference genome
    - **Source:** `genome.txt`
    - **Returns:** Genome string (29,904 nucleotides)

    **2. `parse_mutations(nodeId)`**
    - **Purpose:** Extract mutations for specific phylogenetic node
    - **Source:** `mutations.txt`
    - **Returns:** List of tuples: `(nt_position, original, new, aa_position, aa_change)`

    **3. `construct_variant_genome(genome_seq, mutations)`**
    - **Purpose:** Apply mutations to reference genome
    - **Process:**
      - Convert genome to list
      - For each mutation: Replace nucleotide at position
      - Track applied mutations
    - **Returns:** Variant genome string

    **4. `extract_phylogenetic_features(tree_path)`**
    - **Purpose:** Extract evolutionary features from phylogenetic tree
    - **Source:** `phylogenetic_tree.nwk` (Newick format)
    - **Caching:** Uses HDF5 file `phylo_features_cache.h5`
    - **Features Extracted:**
      - Clade distances (branch lengths)
      - Normalized distances
      - Diversity metrics:
        - Total branch length
        - Max/min branch length
    - **Returns:** Tuple (normalized_distances dict, diversity_metrics dict)

    **5. `calculate_phylogenetic_diversity(tree)`**
    - **Purpose:** Calculate diversity metrics from phylogenetic tree
    - **Metrics:**
      - Total branch length: Sum of all clade depths
      - Max branch length: Deepest clade
      - Min branch length: Shallowest clade
    - **Returns:** Dict of diversity metrics

    **6. `translate_nucleotides_to_amino_acids(nucleotide_sequence, codon_mapper)`**
    - **Purpose:** Convert nucleotide sequence to amino acid sequence
    - **Method:** Sliding window of 3 nucleotides (codons)
    - **Source:** `codon_aa_mapping.json`
    - **Returns:** Amino acid sequence string

    **7. `compute_nucleotide_frequencies(sequence, window_size=10)`**
    - **Purpose:** Calculate local nucleotide composition
    - **Method:** Sliding window frequency calculation
    - **Window:** 10 nucleotides (configurable)
    - **Returns:** List of dicts with frequencies for A, T, G, C

    **8. `preprocess_input(features, expected_size=205, is_multi_input=False)`**
    - **Purpose:** Standardize and encode feature vectors
    - **Process:**
      1. Identify categorical vs. numerical features
      2. One-hot encode categorical features
      3. Standardize numerical features (StandardScaler)
      4. Pad or truncate to expected size (205 features)
      5. Reshape for multi-input models if needed
    - **Returns:** Preprocessed numpy array

    **Precomputation and Caching:**

    **9. `precompute_feature_vectors(...)`**
    - **Purpose:** Generate feature vectors for reference genome
    - **Process for each position:**
      1. Extract k-mer window (30 nucleotides, centered)
      2. Get current nucleotide
      3. Translate to amino acid
      4. Extract phylogenetic features
      5. Compute nucleotide frequencies
      6. Add amino acid biochemical properties:
         - Hydrophobicity
         - Polarity
         - Iso-electric point
         - Volume
         - Molecular weight
      7. Preprocess to 205-dimensional vector
    - **Parameters:**
      - `k=30`: K-mer window size
      - `expected_size=205`: Final feature vector size
    - **Returns:** List of 205-dimensional feature vectors

    **10. `precompute_and_cache_ref_features(...)`**
    - **Purpose:** Cache precomputed features to HDF5
    - **File:** `features.h5`
    - **Cache Structure:**
      - Dataset: "features" (29,904 × 205 matrix)
      - Attributes: genome_length, vector_length
    - **Compression:** GZIP level 9
    - **Behavior:**
      - If cache exists: Load from file
      - If not: Compute and save

    **Variant Processing:**

    **11. `cache_precomputed_features(...)`**
    - **Purpose:** Main feature extraction coordinator
    - **Logic:**
      1. Check cache for reference genome features
      2. If not cached: Call `precompute_and_cache_ref_features()`
      3. If mutations exist:
         - Identify affected positions
         - Call `process_variant()` to update features
      4. Return feature matrix
    - **Cache Key:** Node ID + model parameters
    - **Returns:** (29,904 × 205) feature matrix

    **12. `get_affected_positions(mutations, genome_len, k, protein_regions)`**
    - **Purpose:** Identify positions requiring feature recomputation
    - **Logic:**
      - For each mutation:
        - Calculate k-mer window around mutation
        - Add all positions in window
        - Filter by protein region if specified
    - **Returns:** Sorted set of position indices

    **13. `process_variant(genome_seq, mutations, precomputed_features, ...)`**
    - **Purpose:** Recompute features for mutated regions
    - **Process:**
      1. Load precomputed features
      2. Apply mutations to genome
      3. For each affected position:
         - Recalculate k-mer window
         - Update amino acid translation
         - Recalculate all position-specific features
         - Update feature vector
      4. Add temporal feature (elapsed days)
      5. Add node-specific features (depth, phylogenetic distance)
    - **Optimization:** Only recomputes affected positions
    - **Returns:** Updated feature matrix

    **Utility Functions:**

    **14. `get_aa_features(aa, new_aa, config_file)`**
    - **Purpose:** Extract amino acid biochemical properties
    - **Properties:**
      - Hydrophobicity
      - Polarity
      - Iso-electric point (pI)
      - Volume
      - Molecular weight
      - pKa, pKb, pKx values
      - Isoelectric point (pl)
    - **Returns:** List of 16+ feature values

    **15. `get_sample_depth(depth_file_path, nodeId)`**
    - **Purpose:** Get sequencing depth for phylogenetic node
    - **Source:** `depth_date.json`
    - **Returns:** Depth integer or None

    **16. `balance_dataset(dataset)`**
    - **Purpose:** Balance training data (mutation vs. no mutation)
    - **Method:** Oversampling minority class
    - **Usage:** Model training (not prediction)

    **17. `is_synonymous(current_aa, new_aa)`**
    - **Purpose:** Check if mutation is synonymous (silent)
    - **Returns:** 1 if same amino acid, 0 otherwise

    **18. `find_protein_region(index, protein_regions)`**
    - **Purpose:** Identify which protein region contains a position
    - **Returns:** Protein name or "Non-coding"

---

42. **`genome/configs.py`**
    - **Purpose:** Configuration data for feature extraction
    - **Function:** `configs()` returns dict with:
      - **Protein Regions:** 11 regions with start/end positions
      - **Amino Acid Features:** Dictionaries for each property
        - Polarity: 20 AA values (normalized)
        - Hydrophobicity: 20 AA values
        - Volume: 20 AA values
        - Iso-electric point: 20 AA values
        - Molecular weight: Raw values (0-204)
        - pKa, pKb, pKx: Acid/base properties
        - pI: Isoelectric point
    - **Usage:** Feature engineering and biochemical calculations

43. **`genome/dummyViews.py`**
    - **Purpose:** Testing or development views
    - **Status:** Not used in production

44. **`genome/somecode.py`**
    - **Purpose:** Experimental or utility code
    - **Status:** Likely development artifacts

---

**Data Files**

45. **`genome/genome.txt`**
    - **Purpose:** SARS-CoV-2 reference genome
    - **Format:** FASTA-like (header + sequence)
    - **Length:** 29,904 nucleotides
    - **Source:** Likely Wuhan-Hu-1 reference (NC_045512.2)

46. **`genome/mutations.txt`**
    - **Purpose:** Phylogenetic mutation data
    - **Format:** Tab-separated values
    - **Columns:**
      - Position (1-based)
      - Nucleotide change (REF>ALT)
      - Amino acid position
      - Amino acid change
      - ... (additional metadata)
      - Node ID (column 10)
    - **Usage:** Link phylogenetic nodes to genetic variants

47. **`genome/phylogenetic_tree.nwk`**
    - **Purpose:** Phylogenetic tree of SARS-CoV-2 variants
    - **Format:** Newick format
    - **Content:** Tree structure with branch lengths
    - **Usage:** Extract evolutionary distances and relationships

48. **`genome/codon_aa_mapping.json`**
    - **Purpose:** Genetic code translation table
    - **Format:** JSON dict mapping codons to amino acids
    - **Example:**
      ```json
      {
        "ATG": "M",
        "TAG": "*",
        ...
      }
      ```

49. **`genome/depth_date.json`**
    - **Purpose:** Sequencing depth metadata for phylogenetic nodes
    - **Format:**
      ```json
      {
        "node_id": {
          "depth": 1234,
          "date": "2021-07-03"
        }
      }
      ```

50. **`genome/random_node_ids.txt`**
    - **Purpose:** List of node IDs for testing or sampling
    - **Usage:** Random selection of nodes

---

**Cache and Model Files**

51. **`genome/node_features.h5`** (Generated)
    - **Purpose:** HDF5 cache for precomputed features
    - **Structure:**
      - Groups: Node IDs
      - Datasets: Feature matrices (29,904 × 205)
    - **Benefits:** 
      - Fast loading (no recomputation)
      - Compressed storage

52. **`genome/phylo_features_cache.h5`** (Generated)
    - **Purpose:** HDF5 cache for phylogenetic features
    - **Structure:**
      - Groups: Tree file paths
      - Datasets: 
        - clade_names
        - normalized_distances
        - diversity_metrics

53. **`genome/migrations/`**
    - **Purpose:** Django database migrations
    - **Content:** `__init__.py` (no migrations needed for this app)

---

#### Static Files (`staticfiles/`)

54. **`staticfiles/admin/`**
    - **Purpose:** Django admin interface static files
    - **Content:** CSS, JS, images for admin panel

55. **`staticfiles/rest_framework/`**
    - **Purpose:** Django REST Framework browsable API assets
    - **Content:** CSS, JS, fonts for API browser interface

---

#### Machine Learning Models

56. **`covid19_models/models/` (Referenced but not in tree)**
    - **Location:** `genome_extractor/covid19_models/models/`
    - **Files:** Various `.keras` model files
    - **Examples:**
      - `balanced_data_model.keras` (default)
      - `multi_*_model.keras` (multi-input variants)
    - **Architecture:** Trained TensorFlow/Keras models
    - **Input:** 205-dimensional feature vectors
    - **Output:** 
      - Single-output: Mutation probability (1 value)
      - Dual-output: [P(no mutation), P(mutation)] (2 values)

---

#### Requirements and Dependencies

57. **`requirements.txt`**
    - **Key Dependencies:**
      - Django 5.1.2
      - djangorestframework 3.15.2
      - tensorflow 2.17.0
      - keras 3.6.0
      - numpy 1.26.4
      - pandas 2.2.3
      - scikit-learn 1.5.2
      - biopython 1.84
      - h5py 3.12.1
      - matplotlib (via imports)
      - logomaker (via imports)
    - **Note:** Some packages listed twice (e.g., biopython, Django)

---

### Backend Data Flow

```
API Request
    ↓
views.predict_genome()
    ↓
1. Load Genome & Mutations
    ↓
2. feature_extractor.cache_precomputed_features()
    ├→ Check HDF5 cache
    ├→ If cached: Load features
    └→ If not: Compute features
        ├→ precompute_feature_vectors()
        │   ├→ K-mer extraction
        │   ├→ Phylogenetic features
        │   ├→ Amino acid translation
        │   └→ Biochemical properties
        └→ Cache to HDF5
    ↓
3. Apply Mutations (if any)
    ├→ get_affected_positions()
    └→ process_variant()
        └→ Recompute features for affected positions
    ↓
4. Load ML Model (.keras)
    ↓
5. Run Prediction
    ├→ Single-input: model.predict(features)
    └→ Multi-input: model.predict([features] * 10)
    ↓
6. Post-process Predictions
    ├→ _get_biologically_distributed_probs()
    │   └→ Apply Ti/Tv ratio
    ├→ calculate_genome_data()
    │   └→ Generate probability distributions
    └→ calculate_protein_region_probabilities()
        └→ Aggregate by protein region
    ↓
7. Return JSON Response
    ├→ genomeData: [P(A), P(T), P(G), P(C)] for each position
    ├→ protein_mutation_probs: Average per region
    ├→ genomeSequence: Variant sequence
    └→ metadata
```

---

## 🔄 Full System Data Flow

```
┌─────────────────────────────────────────────────────────────┐
│                        USER ACTIONS                          │
└────────┬────────────────────────────────────────────────────┘
         │
         ├─ Selects Node ID (phylogenetic sample)
         ├─ Enters Elapsed Days (temporal prediction)
         ├─ Selects ML Model
         └─ (Optional) Selects Protein Region
         │
         ↓
┌─────────────────────────────────────────────────────────────┐
│                      FRONTEND (React)                        │
│  Components: Navbar.js → App.js                              │
└────────┬────────────────────────────────────────────────────┘
         │
         ├─ Form Validation
         ├─ Redux State Update (selectNode)
         │
         ↓
         HTTP POST /api/predict/
         {
           nodeId, elapsedDay, selectedModel, selectedProteinRegion
         }
         │
         ↓
┌─────────────────────────────────────────────────────────────┐
│                    BACKEND (Django)                          │
│  views.py → predict_genome()                                 │
└────────┬────────────────────────────────────────────────────┘
         │
         ├─ 1. Load Reference Genome (genome.txt)
         ├─ 2. Parse Mutations (mutations.txt → filter by nodeId)
         ├─ 3. Construct Variant Genome (apply mutations)
         │
         ↓
┌─────────────────────────────────────────────────────────────┐
│              FEATURE EXTRACTION ENGINE                       │
│  feature_extractor.py                                        │
└────────┬────────────────────────────────────────────────────┘
         │
         ├─ Check HDF5 Cache (node_features.h5)
         │   ├─ Cache Hit: Load precomputed features
         │   └─ Cache Miss: Compute features
         │       ├─ K-mer windows (30 nt)
         │       ├─ Phylogenetic distances
         │       ├─ Amino acid translation
         │       ├─ Biochemical properties
         │       └─ Nucleotide frequencies
         │
         ├─ For Variant:
         │   ├─ Identify affected positions
         │   ├─ Recompute features for mutations
         │   └─ Add temporal features (elapsed days)
         │
         └─ Output: Feature Matrix (29,904 × 205)
         │
         ↓
┌─────────────────────────────────────────────────────────────┐
│              MACHINE LEARNING INFERENCE                      │
│  TensorFlow/Keras Model                                      │
└────────┬────────────────────────────────────────────────────┘
         │
         ├─ Load Model: {selectedModel}.keras
         ├─ Input: (29,904 × 205) feature matrix
         ├─ Forward Pass: Neural network inference
         │
         └─ Output: Predictions
             ├─ Single-output: (29,904 × 1) → P(mutation)
             └─ Dual-output: (29,904 × 2) → [P(no_mut), P(mut)]
         │
         ↓
┌─────────────────────────────────────────────────────────────┐
│            PROBABILITY DISTRIBUTION CALCULATION              │
│  views.py → calculate_genome_data()                          │
└────────┬────────────────────────────────────────────────────┘
         │
         ├─ For Each Position:
         │   ├─ Extract prediction values
         │   ├─ Apply biological model (Ti/Tv ratio = 2.0)
         │   │   ├─ Transitions (A↔G, C↔T): Higher probability
         │   │   └─ Transversions: Lower probability
         │   ├─ Normalize probabilities (sum = 1.0)
         │   └─ Output: [P(A), P(T), P(G), P(C)]
         │
         ├─ Aggregate by Protein Region:
         │   ├─ ORF1ab, S, N, M, E, etc.
         │   └─ Calculate average mutation probability
         │
         └─ Output: genomeData + protein_mutation_probs
         │
         ↓
         JSON Response
         {
           genomeData: [[P(A), P(T), P(G), P(C)], ...],  // 29,904 positions
           protein_mutation_probs: {ORF1ab: 0.023, S: 0.045, ...},
           genomeSequence: "ATGC...",
           proteinRegionPossibilities: {ORF1ab: [266, 21555], ...}
         }
         │
         ↓
┌─────────────────────────────────────────────────────────────┐
│                  FRONTEND (React Redux)                      │
│  App.js → setDataset action                                  │
└────────┬────────────────────────────────────────────────────┘
         │
         ├─ Store in Redux (genomeSlice)
         ├─ Navigate to /genome-mutation-visualization
         │
         ↓
┌─────────────────────────────────────────────────────────────┐
│              VISUALIZATION COMPONENTS                        │
└────────┬────────────────────────────────────────────────────┘
         │
         ├─ Recharts.js (Main Chart)
         │   ├─ Render 29,904 bars (or filtered region)
         │   ├─ Color by nucleotide probability
         │   ├─ Add protein region annotations
         │   ├─ Enable zoom/pan
         │   └─ Support WebLogo view for detail
         │
         ├─ DoughnutChart.js (Protein Summary)
         │   ├─ Show mutation % per protein region
         │   └─ Click to filter main chart
         │
         └─ SidePanel.js (Protein Selector)
             ├─ List all protein regions
             └─ Click to zoom to region
         │
         ↓
┌─────────────────────────────────────────────────────────────┐
│                   USER VISUALIZATION                         │
│  - Interactive bar chart                                     │
│  - Protein region highlighting                               │
│  - Mutation hotspot identification                           │
│  - WebLogo for detailed views                                │
└─────────────────────────────────────────────────────────────┘
```

---

## 🧬 Feature Engineering Details

### Feature Vector Composition (205 dimensions)

1. **K-mer Features (30 dimensions)**
   - Sliding window of 30 nucleotides
   - Centered on position of interest
   - One-hot encoded

2. **Position Features (2 dimensions)**
   - Current nucleotide
   - Position index

3. **Phylogenetic Features (1 dimension)**
   - Normalized clade distance from tree

4. **Amino Acid Identity (1 dimension)**
   - Current amino acid at codon position

5. **Nucleotide Frequencies (4 dimensions)**
   - Local composition (A, T, G, C frequencies)
   - Window size: 10 nucleotides

6. **Amino Acid Biochemical Properties (5 dimensions)**
   - Hydrophobicity
   - Polarity
   - Iso-electric point
   - Volume
   - Molecular weight

7. **Additional AA Properties (12 dimensions)**
   - pKa, pKb, pKx values
   - Isoelectric point (pl)

8. **Mutation-Specific Features (when applicable)**
   - Reference vs. alternate amino acid
   - Synonymous/non-synonymous indicator
   - Protein region indicator

9. **Temporal Features**
   - Elapsed days (for temporal prediction)

10. **Padding** (to reach 205 dimensions)
    - Zero-padding for consistent vector size

---

## 🧠 Machine Learning Models

### Model Types

1. **Single-Input Models**
   - Input: (batch_size, 205)
   - Output: (batch_size, 1) or (batch_size, 2)
   - Example: `balanced_data_model.keras`

2. **Multi-Input Models**
   - Input: (batch_size, 10, 205)
   - Output: (batch_size, 1) or (batch_size, 2)
   - Architecture: Likely LSTM or Conv1D
   - Example: Models with "multi" in filename

### Output Interpretation

- **Single-Output (shape: [n, 1])**
  - Direct mutation probability: P(mutation)
  - P(no mutation) = 1 - P(mutation)

- **Dual-Output (shape: [n, 2])**
  - Explicit: [P(no mutation), P(mutation)]
  - More robust probabilistic outputs

### Biological Post-Processing

**Transition/Transversion (Ti/Tv) Model:**
- Purines: A, G
- Pyrimidines: T, C
- Transitions (A↔G or T↔C): More common (2× likely)
- Transversions (purine↔pyrimidine): Less common

Formula:
```
prob_tv = P(mutation) / (ti_tv_ratio + 2)
prob_ti = prob_tv * ti_tv_ratio
```

Where `ti_tv_ratio = 2.0`

---

## 🔧 Key Technologies and Libraries

### Frontend Libraries

| Library | Version | Purpose |
|---------|---------|---------|
| React | 18.2.0 | UI framework |
| Redux Toolkit | 1.9.7 | State management |
| React Router | 6.14.2 | Routing |
| Chart.js | 4.4.0 | Charting engine |
| chartjs-plugin-zoom | 2.0.1 | Chart zoom/pan |
| chartjs-plugin-annotation | 3.0.1 | Chart annotations |
| Recharts | 2.10.3 | Alternative charting |
| Material Tailwind | 2.1.10 | UI components |
| Tailwind CSS | 4.1.7 | Styling framework |
| Axios | 1.6.3 | HTTP client |
| D3.js | 7.8.5 | Data visualization |
| Plotly.js | 2.27.1 | Scientific plotting |

### Backend Libraries

| Library | Version | Purpose |
|---------|---------|---------|
| Django | 5.1.2 | Web framework |
| DRF | 3.15.2 | REST API |
| TensorFlow | 2.17.0 | ML inference |
| Keras | 3.6.0 | Model interface |
| NumPy | 1.26.4 | Numerical computing |
| Pandas | 2.2.3 | Data manipulation |
| Scikit-Learn | 1.5.2 | Preprocessing |
| Biopython | 1.84 | Biological data |
| h5py | 3.12.1 | HDF5 caching |
| Matplotlib | (latest) | Plotting |
| Logomaker | (latest) | Sequence logos |

---

## 📊 Performance Optimizations

### Frontend
1. **Redux State Management:** Centralized, prevents unnecessary re-renders
2. **Virtualized Dropdowns:** React-select-virtualized for large lists
3. **Canvas Rendering:** Chart.js uses canvas (faster than SVG)
4. **Conditional Rendering:** High-res views only when needed
5. **Memoization:** useMemo hooks for expensive calculations
6. **Code Splitting:** React Router lazy loading (if implemented)

### Backend
1. **HDF5 Caching:** Precomputed features stored in compressed HDF5
2. **Selective Recomputation:** Only mutated regions reprocessed
3. **Batch Prediction:** All positions predicted in single forward pass
4. **NumPy Vectorization:** Avoid Python loops where possible
5. **Django Caching:** Session-based caching for genome data
6. **GZIP Compression:** HDF5 datasets compressed (level 9)

---

## 🔒 Security Considerations

### Current Issues
1. **Secret Key Exposed:** In `settings.py` (should be in env vars)
2. **Debug Mode Enabled:** `DEBUG = True` in settings
3. **CORS Wide Open:** Allows localhost, but needs production review
4. **No Authentication:** API endpoints are open (consider adding auth)
5. **SQL Injection:** Minimal risk (using Django ORM, but models.py empty)

### Recommendations
1. Move secret key to environment variables
2. Disable debug mode in production
3. Implement API authentication (JWT, API keys)
4. Add rate limiting for API endpoints
5. Validate and sanitize all user inputs
6. Use HTTPS in production
7. Implement CSRF protection for POST requests

---

## 🚀 Deployment Notes

### Frontend Deployment
- **Build Command:** `npm run build`
- **Output:** `build/` directory with static files
- **Hosting Options:** 
  - Nginx (static files)
  - Vercel, Netlify (automatic deployment)
  - AWS S3 + CloudFront

### Backend Deployment
- **WSGI Server:** Gunicorn, uWSGI
- **Web Server:** Nginx (reverse proxy)
- **Database:** SQLite (dev), PostgreSQL (production recommended)
- **Static Files:** `python manage.py collectstatic`
- **Domain:** Currently configured for `covidmutex.itu.edu.tr`

### Production Checklist
- [ ] Set `DEBUG = False`
- [ ] Use environment variables for secrets
- [ ] Configure production database
- [ ] Set up HTTPS/SSL
- [ ] Configure CORS properly
- [ ] Set up logging and monitoring
- [ ] Implement backup strategy for HDF5 caches
- [ ] Optimize model loading (keep in memory)
- [ ] Set up CDN for static files
- [ ] Implement health check endpoints

---

## 📈 Scalability Considerations

### Current Bottlenecks
1. **Model Loading:** Loads .keras file on each request
2. **Feature Extraction:** CPU-intensive for uncached nodes
3. **Single-threaded Django:** Default development server
4. **In-memory Processing:** 29,904-position arrays

### Scalability Solutions
1. **Model Caching:** Load models once, keep in memory
2. **Worker Processes:** Gunicorn with multiple workers
3. **Async Processing:** Celery for long-running tasks
4. **Database:** Move from SQLite to PostgreSQL
5. **Redis Caching:** Cache API responses
6. **GPU Inference:** Use TensorFlow GPU for faster predictions
7. **CDN:** Cache static assets and frequently-used API responses
8. **Load Balancer:** Distribute traffic across multiple servers

---

## 🧪 Testing Strategy

### Frontend Testing (Recommended)
- **Unit Tests:** Jest + React Testing Library
- **Component Tests:** Test individual components
- **Integration Tests:** Test Redux integration
- **E2E Tests:** Cypress or Playwright

### Backend Testing (Recommended)
- **Unit Tests:** Django TestCase
- **API Tests:** DRF APITestCase
- **Feature Extraction Tests:** Validate feature vectors
- **Model Tests:** Validate model outputs
- **Performance Tests:** Benchmark feature extraction

---

## 📚 Domain Knowledge

### SARS-CoV-2 Genome Structure
- **Length:** 29,904 nucleotides
- **Type:** Positive-sense single-stranded RNA
- **Protein Regions:**
  - ORF1ab: Replicase polyprotein (largest, 21,290 nt)
  - S: Spike protein (3,822 nt) - vaccine target
  - N: Nucleocapsid protein (1,260 nt)
  - M: Membrane protein (669 nt)
  - E: Envelope protein (228 nt)
  - ORF3a, ORF6, ORF7a, ORF7b, ORF8, ORF10: Accessory proteins

### Mutation Types
- **Synonymous:** Silent mutation, same amino acid
- **Non-synonymous:** Changes amino acid
- **Transition:** Purine↔Purine or Pyrimidine↔Pyrimidine
- **Transversion:** Purine↔Pyrimidine
- **Indel:** Insertion or deletion (not modeled here)

### Phylogenetic Analysis
- **Node:** Represents a sampled viral sequence
- **Branch Length:** Evolutionary distance
- **Clade:** Group of related samples
- **Tree:** Represents evolutionary relationships

---

## 🎯 Use Cases

1. **Vaccine Development:**
   - Identify mutation hotspots in Spike protein
   - Predict future variants
   - Design broad-spectrum vaccines

2. **Therapeutic Target Identification:**
   - Find conserved regions (low mutation probability)
   - Identify stable drug targets

3. **Evolutionary Analysis:**
   - Track mutation patterns over time
   - Understand selective pressures

4. **Public Health:**
   - Monitor emerging variants
   - Assess variant impact on transmissibility

---

## 🐛 Known Issues and TODOs

### Known Issues
1. No package.json in frontend root (only .backup)
2. Duplicate dependencies in requirements.txt
3. Empty models.py (database not utilized)
4. Multiple chart libraries (could consolidate)
5. No error logging/monitoring
6. No API rate limiting
7. No user authentication

### TODOs
1. Add comprehensive testing
2. Implement API authentication
3. Add API documentation (Swagger/OpenAPI)
4. Optimize model loading
5. Add logging and monitoring
6. Implement backup strategy
7. Add user accounts and saved analyses
8. Export functionality (CSV, JSON, images)
9. Batch analysis support
10. Comparison between multiple nodes

---

## 📖 Documentation

### Files Needing Documentation
- API endpoint documentation
- Feature extraction algorithm details
- Model training methodology
- Data preprocessing pipeline
- Deployment instructions
- User manual

### Recommended Tools
- **API Docs:** drf-spectacular (OpenAPI/Swagger)
- **Code Docs:** JSDoc (frontend), Sphinx (backend)
- **Architecture Diagrams:** Draw.io, PlantUML
- **User Guide:** Markdown with screenshots

---

## 🤝 Contributing

### Development Workflow
1. **Setup:**
   - Frontend: `cd covid19-genome-visualizer && npm install`
   - Backend: `cd genome_extractor && pip install -r requirements.txt`

2. **Run Development Servers:**
   - Frontend: `npm start` (port 3000)
   - Backend: `python manage.py runserver` (port 8000)

3. **Make Changes:**
   - Follow existing code style
   - Add tests for new features
   - Update documentation

4. **Commit:**
   - Use descriptive commit messages
   - Reference issues if applicable

### Code Style
- **Frontend:** Prettier, ESLint
- **Backend:** Black, Flake8, isort

---

## 📊 Project Statistics

- **Total Lines of Code:** ~15,000+ (estimated)
- **Frontend Components:** 27
- **Backend Views:** 2 main endpoints
- **Feature Dimensions:** 205
- **Genome Length:** 29,904 nucleotides
- **Protein Regions:** 11
- **ML Models:** Multiple (.keras files)
- **Dependencies:** 40+ libraries

---

## 🏆 Key Achievements

1. **Integration of ML and Bioinformatics:** Combines deep learning with genomic analysis
2. **Interactive Visualization:** Multiple chart types for comprehensive data exploration
3. **Performance:** Efficient caching and vectorization
4. **Biological Accuracy:** Ti/Tv ratio modeling for realistic predictions
5. **Scalable Architecture:** Modular design allows for easy extension
6. **User-Friendly Interface:** Intuitive controls and responsive design

---

## 🔮 Future Enhancements

1. **Real-time Updates:** WebSocket for live predictions
2. **Comparative Analysis:** Side-by-side comparison of multiple variants
3. **Export Functionality:** PDF reports, CSV data export
4. **Advanced Filtering:** Filter by mutation type, impact, protein function
5. **Mobile Support:** Responsive design optimization
6. **Collaboration Features:** Share analyses, annotations
7. **Integration with Public Databases:** GISAID, NCBI, etc.
8. **Automated Model Retraining:** Periodic updates with new data
9. **Multi-language Support:** Internationalization
10. **Advanced Visualizations:** 3D protein structure, phylogenetic tree viewer

---

## 📝 Conclusion

CovMutEx is a sophisticated bioinformatics tool that successfully combines machine learning, genomic analysis, and interactive visualization to address a critical need in pandemic research. The architecture is well-designed with clear separation of concerns, though there are opportunities for optimization and enhancement, particularly around security, testing, and scalability.

The project demonstrates strong technical competence across multiple domains:
- Full-stack web development (React + Django)
- Machine learning (TensorFlow/Keras)
- Bioinformatics (feature extraction, phylogenetics)
- Data visualization (multiple charting libraries)
- Scientific computing (NumPy, Pandas, BioPython)

With some refinements and additional features, this tool has the potential to become a valuable resource for the virology and vaccine development community.

---

**Last Updated:** October 8, 2025
**Repository:** itu-bioinformatics-database-lab/CovMutEx
**Branch:** clean-main
