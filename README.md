MATLAB implementations for kidney and lung kinetic modeling as local add-ons within the Turku PET Centre **TURBO** pipeline.  
These tools provide automated ROI extraction, mask generation, TAC extraction, and kinetic parameter estimation for **[15O]H2O** PET imaging.

The repository currently includes:

- **`turbo_kidney.m`** — Kidney kinetic modeling
- **`turbo_lung.m`** — Lung kinetic modeling
- **`turbo_kidney_petmr.m`** — Kidney kinetic modeling using pet/mr instead of Quadra
- **`run_rhIDIF.m`** — Script for using RH IDIF for extracting input function (https://github.com/Rigshospitalet-KFNM/IDIF)
- **`orthanc_export_worker.py`** — Python script for moving incomming DICOM files to turbo file structure using Orthanc dicom node and LUA
 
All scripts follow the same directory structure as the official **TURBO** (Turku Image-Based Input Function Optimization) framework (https://turbo.utu.fi/)

------------------------------------------------------------------------
# Orthanc Export Worker

Automated DICOM Export Pipeline (Orthanc + Lua + Python)

## Overview

This project implements an automated DICOM export pipeline using: -
Orthanc (DICOM server) - Lua hook (event trigger) - Python export worker
(background daemon)

When Orthanc receives DICOM instances, a Lua script creates lightweight
job files.\
A Python worker monitors these jobs, waits until series are stable, and
exports DICOM files into a structured folder hierarchy.

Key features: - Non-blocking ingestion (Orthanc stays fast) - Robust for
very large series (10k--30k+ images) - Idempotent (no duplicate
exports) - Crash-safe and restart-safe

------------------------------------------------------------------------

## System Architecture

    DICOM Sender (Scanner / PACS)
            │
            ▼
          Orthanc (C-STORE)
            │
            ▼
       Lua OnStoredInstance Hook
            │
            ▼
     Job Queue (/var/lib/orthanc/jobs)
            │
            ▼
     Python Export Worker (systemd service)
            │
            ▼
     Structured Export Folder
     /home/turbo/Documents/data/<ProjectID>/dicom/...

------------------------------------------------------------------------

## Components

### 1. Orthanc (DICOM Server)

Orthanc receives and stores DICOM data and exposes a REST API used by
the worker.

Important directories: - Storage: `/var/lib/orthanc/db-v6` - Jobs:
`/var/lib/orthanc/jobs`

Orthanc REST endpoints used: - `/instances/{id}` - `/series/{id}` -
`/studies/{id}` - `/series/{id}/archive`

------------------------------------------------------------------------

### 2. Lua Hook (on_receive.lua)

Location:

    /etc/orthanc/lua/on_receive.lua

#### Purpose

Triggered on every stored instance (`OnStoredInstance`).\
Instead of heavy processing inside Orthanc, it writes a small job file:

    /var/lib/orthanc/jobs/<instanceId>.job

Example content:

    instanceId=<OrthancInstanceUUID>

Why this design: - Prevents blocking Orthanc - Avoids crashes during
large exports - Enables asynchronous processing

------------------------------------------------------------------------

### 3. Python Export Worker

Service:

    orthanc-export-worker.service

Script:

    /usr/local/bin/orthanc_export_worker.py

The worker: - Watches the job directory - Maps instance → series →
study - Tracks stability of incoming data - Exports DICOM files when
series are complete

------------------------------------------------------------------------

## Processing Workflow (Step-by-Step)

### Step 1 --- Job Detection

The worker scans:

    /var/lib/orthanc/jobs/*.job

It: 1. Atomically renames `.job` → `.working` 2. Reads `instanceId` 3.
Deletes the working file

This prevents race conditions and duplicate processing.

------------------------------------------------------------------------

### Step 2 --- Instance Mapping

Using Orthanc REST:

    GET /instances/{instanceId}

The worker retrieves: - ParentSeries - ParentStudy

Fallback logic ensures compatibility across Orthanc versions.

------------------------------------------------------------------------

### Step 3 --- State Tracking

State file:

    /var/lib/orthanc/jobs/study_state.json

Stored data: - Last time a study received data - Per-series instance
counts - Stability timestamps

Purpose: - Avoid exporting incomplete series - Support delayed
multi-series studies - Enable safe restarts

------------------------------------------------------------------------

### Step 4 --- Series Stability Detection

Large series (e.g., 30,000 images) may arrive over minutes or hours.

Instead of a fixed timeout, the worker: 1. Queries `/series/{id}` 2.
Reads `NumberOfInstances` 3. Waits until the count stops increasing
for: - `SERIES_QUIET_SECONDS` (default: 30s) - OR
`SERIES_MAX_WAIT_SECONDS` (safety cap)

This prevents partial exports.

------------------------------------------------------------------------

### Step 5 --- Memory-Safe Export (Streaming)

Export endpoint:

    GET /series/{seriesId}/archive

The worker: - Streams the ZIP archive to a temporary file (NOT RAM) -
Extracts DICOM files - Deletes the temporary ZIP

Temporary data locations: - Temp ZIP: same filesystem as destination (or
configurable) - Final DICOM: destination folder only

This avoids Out-Of-Memory crashes with large series.

------------------------------------------------------------------------

### Step 6 --- Flattened Extraction

Orthanc archives contain nested folders:

    Patient/Study/Series/*.dcm

The worker flattens extraction:

    CT/
      IMG0001.dcm
      IMG0002.dcm

No nested directories are kept.

------------------------------------------------------------------------

## Export Folder Structure

### Base Path

    /home/turbo/Documents/data/<ProjectID>/dicom/<PatientID>_<StudyDate>_<state>/<MODALITY>

### Examples

PET:

    .../448/dicom/12345_20260216_vaagen/PET

CT:

    .../448/dicom/12345_20260216_vaagen/CT

Special rule: - If PET SeriesDescription contains `" O2 "` → `_o2` is
added:

    .../12345_20260216_o2_vaagen/PET

State is automatically lowercased.

------------------------------------------------------------------------

## Study Description Parsing

Example:

    Specials^448_IMPACT_Pressor_1 (Adult)

Parsed as: - ProjectID = 448 - state = pressor1

The parser is robust to inconsistent naming and extra underscores.

------------------------------------------------------------------------

## Idempotency (No Duplicate Exports)

Each exported series creates a marker file:

    .exported_<SeriesOrthancID>.done

Benefits: - Safe restarts - No duplicate downloads - Retry-safe
processing

------------------------------------------------------------------------

## Temporary Data Summary

  ------------------------------------------------------------------------------------------
  Component                Location                                   Purpose
  ------------------------ ------------------------------------------ ----------------------
  Orthanc storage          `/var/lib/orthanc/db-v6`                   Raw DICOM database

  Job queue                `/var/lib/orthanc/jobs`                    Trigger files

  State file               `/var/lib/orthanc/jobs/study_state.json`   Processing state

  Temp ZIP (export)        Destination filesystem                     Streaming archive

  Final output             `/home/turbo/Documents/data`               Exported DICOM
  ------------------------------------------------------------------------------------------

------------------------------------------------------------------------

## Performance Considerations (Large Series)

Recommended: - SSD storage - ≥16 GB RAM - Sufficient disk space for
export + temporary ZIP files

Critical design choices: - Streaming download (prevents OOM) -
Per-series stability detection - Asynchronous worker (non-blocking
Orthanc)

------------------------------------------------------------------------

## Running the Services

Start Orthanc:

    sudo systemctl start orthanc

Start export worker:

    sudo systemctl start orthanc-export-worker

Check logs:

    journalctl -u orthanc-export-worker -n 100 --no-pager



------------------------------------------------------------------------

## Design Rationale

Heavy export tasks are intentionally NOT performed inside Lua/Orthanc
because: - Lua runs inside the Orthanc process - Large exports would
block DICOM ingestion - Increased crash risk

The job-queue + worker architecture provides: - Isolation -
Reliability - Scalability - Clinical-grade robustness for large datasets
