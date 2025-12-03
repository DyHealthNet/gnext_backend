# GNExT Platform - Backend

<div align="center">
  <img src="src/assets/figures/GNExT_Logo_Black.png" alt="GNExT Logo" width="400"/>
  <p><em>GWAS Network Exploration Tool</em></p>
</div>

## Overview

The backend of the GNExT platform is a Django REST API providing data access, search capabilities, and processing services for GWAS (Genome-Wide Association Studies) data analysis and visualization.

For detailed overview of the platform, please refer to the **[GNExT platform repository](https://github.com/DyHealthNet/gnext_platform)** which describes the platform in more detail and showcases how to set up a GNExT platform for your study data. The repository provides information on the configuration parameters and the deployment strategy.

### 🔗 Repository Links

- **📦 [GNExT Platform](https://github.com/DyHealthNet/gnext_platform)** - Complete platform with frontend, backend, and deployment
- **🔧 [GNExT Backend](https://github.com/DyHealthNet/gnext_backend)** - Django REST API and data processing (this repository)
- **🎨 [GNExT Frontend](https://github.com/DyHealthNet/gnext_frontend)** - Vue.js web interface
- **⚙️ [GNExT Nextflow Pipeline](https://github.com/DyHealthNet/gnext_nf_pipeline)** - Data processing pipeline for GWAS analysis


## Development Setup

In development mode, ensure that the [GNExT Platform](https://github.com/DyHealthNet/gnext_platform) repository is cloned with all submodules. After completing the .env configuration, proceed to the backend directory to initialize and configure the development environment according to the steps outlined below.

### Quick Start

**Option 1: Create new conda environment**
```bash
# Create and activate conda environment with Python 3.11
conda create -n gnext_backend python=3.11
conda activate gnext_backend

# Install Python dependencies
pip install -r requirements.txt

# Initialize Typesense (if needed)
python init_typesense.py

# Start development server
python manage.py runserver 0.0.0.0:8300
```

**Option 2: Use existing conda environment**
```bash
# Activate your existing conda environment
conda activate gnext_backend

# Install/update Python dependencies
pip install -r requirements.txt

# Initialize Typesense (if needed)
python init_typesense.py

# Start development server
python manage.py runserver 0.0.0.0:8300
```

The backend API will be available at `http://localhost:8300` (or the port specified in `VITE_BACKEND_PORT`)

## Deployment Setup

Please have a look at the [GNExT Platform](https://github.com/DyHealthNet/gnext_platform) repository.


## Citation
TBA