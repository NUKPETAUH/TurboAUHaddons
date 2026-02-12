MATLAB implementations for kidney and lung kinetic modeling as local add-ons within the Turku PET Centre **TURBO** pipeline.  
These tools provide automated ROI extraction, mask generation, TAC extraction, and kinetic parameter estimation for **[15O]H2O** PET imaging.

The repository currently includes:

- **`turbo_kidney.m`** — Kidney kinetic modeling
- **`turbo_lung.m`** — Lung kinetic modeling
- **`turbo_kidney_petmr.m`** — Kidney kinetic modeling using pet/mr instead of Quadra
- **`run_rhIDIF.m`** — Script for using RH IDIF for extracting input function (https://github.com/Rigshospitalet-KFNM/IDIF)
- **`orthanc_export_worker.py`** — Python script for moving incomming DICOM files to turbo file structure using Orthanc dicom node and LUA
 
All scripts follow the same directory structure as the official **TURBO** (Turku Image-Based Input Function Optimization) framework (https://turbo.utu.fi/)
