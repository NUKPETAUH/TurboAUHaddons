#!/usr/bin/env python3
import os, time, json, re, zipfile, io
import requests
import tempfile
import shutil

# --------- CONFIG ----------
ORTHANC_URL = "http://127.0.0.1:8042"
ORTHANC_USER = ""
ORTHANC_PASS = ""

JOB_DIR = "/var/lib/orthanc/jobs"
STATE_FILE = os.path.join(JOB_DIR, "study_state.json")

# Tuning for large series
QUIET_SECONDS = 20                  # small study-level quiet
SERIES_QUIET_SECONDS = 30           # series count must be stable for this long
SERIES_MAX_WAIT_SECONDS = 2 * 60 * 60  # safety cap (2 hours)

BASE_OUT = "/home/turbo/Documents/data"
# --------------------------


def make_session():
    s = requests.Session()
    s.headers.update({"Connection": "close"})  # Orthanc keep-alive can be short
    if ORTHANC_USER:
        s.auth = (ORTHANC_USER, ORTHANC_PASS)
    return s


sess = make_session()


def rest_get(path: str):
    global sess
    try:
        r = sess.get(ORTHANC_URL + path, timeout=30)
        r.raise_for_status()
        return r.json()
    except (requests.exceptions.ConnectionError, requests.exceptions.ChunkedEncodingError):
        sess = make_session()
        r = sess.get(ORTHANC_URL + path, timeout=30)
        r.raise_for_status()
        return r.json()


def rest_get_bytes(path: str):
    global sess
    try:
        r = sess.get(ORTHANC_URL + path, timeout=300)
        r.raise_for_status()
        return r.content
    except (requests.exceptions.ConnectionError, requests.exceptions.ChunkedEncodingError):
        sess = make_session()
        r = sess.get(ORTHANC_URL + path, timeout=300)
        r.raise_for_status()
        return r.content


def parse_study_description(sd: str):

    """
    Handles inconsistent formats like:
      Specials^448_IMPACT_Vaagen (Adult)
      Specials^448_IMPACT_Pressor_1 (Adult)
      Specials^448_SOMETHING_Long_State_Name (Adult)

    Returns:
      ProjectID = 448
      state = Pressor1 / Vaagen / LongStateName
    """
    if not sd:
        return ("UNK", "UNK")

    # Extract ProjectID (digits after ^)
    pid_match = re.search(r"\^(\d+)", sd)
    project_id = pid_match.group(1) if pid_match else "UNK"

    # Remove the prefix up to first caret and any trailing parentheses
    # Example: "Specials^448_IMPACT_Pressor_1 (Adult)" → "448_IMPACT_Pressor_1"
    core = re.sub(r"\s*\(.*\)$", "", sd)  # remove "(Adult)"
    if "^" in core:
        core = core.split("^", 1)[1]

    # Split by underscore
    parts = core.split("_")

    # Expected pattern: [ProjectID, SOMETHING, STATE...]
    if len(parts) >= 3:
        state_parts = parts[2:]
        state_raw = "".join(state_parts)  # join Pressor + 1 → Pressor1
    else:
        state_raw = "UNK"

    # Keep alphanumeric only (safe for folders)
    state = re.sub(r"[^A-Za-z0-9]+", "", state_raw)

    return (project_id, state if state else "UNK")

def safe(s: str):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s or "UNK")


def contains_o2(series_desc: str) -> bool:
    if not series_desc:
        return False
    # match your exact rule: contains ' O2 ' (with spaces), case-insensitive
    return " O2 " in series_desc.upper()


def ensure_dir(p: str):
    os.makedirs(p, exist_ok=True)


def load_state():
    if os.path.exists(STATE_FILE):
        with open(STATE_FILE, "r") as f:
            return json.load(f)
    return {}


def save_state(state):
    tmp = STATE_FILE + ".tmp"
    with open(tmp, "w") as f:
        json.dump(state, f)
    os.replace(tmp, STATE_FILE)


def series_marker_path(dest_dir: str, series_id: str) -> str:
    return os.path.join(dest_dir, f".exported_{series_id}.done")


def is_series_exported(dest_dir: str, series_id: str) -> bool:
    return os.path.exists(series_marker_path(dest_dir, series_id))


def mark_series_exported(dest_dir: str, series_id: str):
    ensure_dir(dest_dir)
    tmp = series_marker_path(dest_dir, series_id) + ".tmp"
    final = series_marker_path(dest_dir, series_id)
    with open(tmp, "w") as f:
        f.write("ok\n")
    os.replace(tmp, final)



def export_series_to(series_id: str, dest_dir: str):
    """
    Stream Orthanc /series/{id}/archive to a temp zip file on disk,
    then extract files into dest_dir (flattened), then delete temp.
    """
    ensure_dir(dest_dir)

    # Put temp zip on the SAME filesystem as dest_dir to avoid cross-device issues
    # (or use /var/tmp if you prefer)
    tmp_dir = dest_dir
    fd, tmp_zip = tempfile.mkstemp(prefix=f"series_{series_id}_", suffix=".zip", dir=tmp_dir)
    os.close(fd)

    try:
        # Stream download to disk (no .content)
        with sess.get(ORTHANC_URL + f"/series/{series_id}/archive", stream=True, timeout=300) as r:
            r.raise_for_status()
            with open(tmp_zip, "wb") as f:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        f.write(chunk)

        # Extract flattened
        with zipfile.ZipFile(tmp_zip) as z:
            for info in z.infolist():
                if info.is_dir():
                    continue

                filename = os.path.basename(info.filename)
                if not filename:
                    continue

                out_path = os.path.join(dest_dir, filename)

                if os.path.exists(out_path):
                    base, ext = os.path.splitext(filename)
                    i = 1
                    while True:
                        candidate = os.path.join(dest_dir, f"{base}_{i}{ext}")
                        if not os.path.exists(candidate):
                            out_path = candidate
                            break
                        i += 1

                with z.open(info) as src, open(out_path, "wb") as dst:
                    shutil.copyfileobj(src, dst, length=1024 * 1024)

    finally:
        try:
            os.remove(tmp_zip)
        except FileNotFoundError:
            pass



def study_metadata(study_id: str):
    """
    Return (patient_id, study_date, study_description).
    PatientID may not be present at study-level; fallback to ParentPatient.
    """
    s = rest_get(f"/studies/{study_id}")
    study_tags = (s.get("MainDicomTags") or {})

    patient_id = study_tags.get("PatientID") or ""
    if not patient_id:
        patient_ref = s.get("ParentPatient")
        if patient_ref:
            p = rest_get(f"/patients/{patient_ref}")
            p_tags = (p.get("MainDicomTags") or {})
            patient_id = p_tags.get("PatientID") or ""

    if not patient_id:
        patient_id = "UNKNOWN"

    study_date = study_tags.get("StudyDate", "UNKNOWN")
    study_desc = study_tags.get("StudyDescription", "")
    return patient_id, study_date, study_desc


def instance_to_series_and_study(instance_id: str):
    """
    Return (series_id, study_id) for an instance, robustly.
    """
    inst = rest_get(f"/instances/{instance_id}")
    series_id = inst.get("ParentSeries")
    study_id = inst.get("ParentStudy")

    if not study_id and series_id:
        series = rest_get(f"/series/{series_id}")
        study_id = series.get("ParentStudy")

    return series_id, study_id


def series_info(series_id: str):
    """
    Return (modality, series_description, number_of_instances)
    """
    s = rest_get(f"/series/{series_id}")
    tags = (s.get("MainDicomTags") or {})
    modality = (tags.get("Modality") or "NA").upper()
    series_desc = tags.get("SeriesDescription") or ""
    num = int(s.get("NumberOfInstances") or 0)
    return modality, series_desc, num


def compute_dest(study_id: str, series_id: str):
    """
    Compute destination directory for this series based on rules.
    Returns (dest_dir, modality_subfolder) where dest_dir is full path to subfolder.
    """
    patient_id, study_date, study_desc = study_metadata(study_id)
    project_id, state = parse_study_description(study_desc)
    state_lower = (state or "UNK").lower()

    modality, series_desc, _ = series_info(series_id)

    if modality in ("PT", "PET"):
        sub = "PET"
        o2_suffix = "_o2" if contains_o2(series_desc) else ""
    elif modality == "CT":
        sub = "CT"
        o2_suffix = ""
    elif modality in ("MR", "MRI"):
        sub = "MRI"
        o2_suffix = ""
    else:
        return None, None  # unsupported modality

    folder_name = f"{safe(patient_id)}_{safe(study_date)}{o2_suffix}_{safe(state_lower)}"
    root = os.path.join(BASE_OUT, safe(project_id), "dicom", folder_name)
    dest = os.path.join(root, sub)
    return dest, sub


def claim_job(job_path: str):
    working = job_path[:-4] + ".working"
    try:
        os.replace(job_path, working)  # atomic
        return working
    except FileNotFoundError:
        return None


def read_instance_id(jobfile: str):
    with open(jobfile, "r") as f:
        for line in f:
            if line.startswith("instanceId="):
                return line.strip().split("=", 1)[1].strip()
    return None


def main():
    ensure_dir(JOB_DIR)
    state = load_state()

    while True:
        now = time.time()

        # 1) consume jobs
        for fn in os.listdir(JOB_DIR):
            if not fn.endswith(".job"):
                continue

            claimed = claim_job(os.path.join(JOB_DIR, fn))
            if not claimed:
                continue

            instance_id = read_instance_id(claimed)
            try:
                os.unlink(claimed)
            except FileNotFoundError:
                pass

            if not instance_id:
                continue

            try:
                series_id, study_id = instance_to_series_and_study(instance_id)
            except Exception as e:
                print(f"worker: ERROR mapping instance {instance_id}: {e}", flush=True)
                continue

            if not series_id or not study_id:
                print(f"worker: could not map instance {instance_id} -> series/study (series={series_id}, study={study_id})", flush=True)
                continue

            rec = state.get(study_id, {"last_seen": 0, "series": {}})
            rec["last_seen"] = now

            # update per-series tracking
            series_rec = rec["series"].get(series_id)
            try:
                _, _, current_count = series_info(series_id)
            except Exception:
                current_count = series_rec.get("last_count", 0) if series_rec else 0

            if series_rec is None:
                rec["series"][series_id] = {
                    "first_seen": now,
                    "last_seen": now,
                    "last_count": current_count,
                    "last_count_changed": now
                }
            else:
                if current_count != series_rec.get("last_count", 0):
                    series_rec["last_count"] = current_count
                    series_rec["last_count_changed"] = now
                series_rec["last_seen"] = now
                rec["series"][series_id] = series_rec

            state[study_id] = rec

        # 2) export series that look stable
        changed = False

        for study_id, rec in list(state.items()):
            if now - rec.get("last_seen", 0) < QUIET_SECONDS:
                continue

            # iterate over tracked series for this study
            for series_id, series_rec in list(rec.get("series", {}).items()):
                dest, sub = compute_dest(study_id, series_id)
                if dest is None:
                    # unsupported modality; drop it from tracking
                    rec["series"].pop(series_id, None)
                    changed = True
                    continue

                if is_series_exported(dest, series_id):
                    rec["series"].pop(series_id, None)
                    changed = True
                    continue

                # refresh instance count
                try:
                    _, _, latest_count = series_info(series_id)
                except Exception as e:
                    print(f"worker: WARNING could not read series {series_id} count: {e}", flush=True)
                    latest_count = series_rec.get("last_count", 0)

                # update stability record if count changed
                if latest_count != series_rec.get("last_count", 0):
                    series_rec["last_count"] = latest_count
                    series_rec["last_count_changed"] = now
                    series_rec["last_seen"] = now
                    rec["series"][series_id] = series_rec
                    changed = True
                    continue

                stable_seconds = now - series_rec.get("last_count_changed", now)
                age_seconds = now - series_rec.get("first_seen", now)

                if stable_seconds >= SERIES_QUIET_SECONDS or age_seconds >= SERIES_MAX_WAIT_SECONDS:
                    print(f"worker: exporting series {series_id} -> {dest} (stable={stable_seconds:.0f}s age={age_seconds:.0f}s count={latest_count})", flush=True)
                    try:
                        export_series_to(series_id, dest)
                        mark_series_exported(dest, series_id)
                        rec["series"].pop(series_id, None)
                        changed = True
                    except Exception as e:
                        print(f"worker: ERROR exporting series {series_id}: {e}", flush=True)
                        # keep it for retry

            state[study_id] = rec

            # optional: drop empty studies to keep state small
            if not state[study_id].get("series"):
                # keep last_seen but remove study if you want; here we remove it
                state.pop(study_id, None)
                changed = True

        if changed:
            save_state(state)

        time.sleep(1)


if __name__ == "__main__":
    main()
