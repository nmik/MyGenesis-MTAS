#!/usr/bin/env python3
"""Local web UI: upload Qβ + γ XML, inspect, and plot the level ladder."""

from __future__ import annotations

import base64
import io
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np
from flask import Flask, jsonify, request, send_from_directory

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from plotLevelLadder import plot_ladder  # noqa: E402
from xmlToLevelLadder import parse_nuclide_xml, run_translation  # noqa: E402

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 16 * 1024 * 1024

SAMPLE_QBETA = ROOT / "98ZR.xml"
SAMPLE_GAMMA = ROOT / "98NB.xml"


def _write_upload(storage, suffix=".xml"):
    raw = storage.read()
    name = Path(storage.filename or "upload.xml").name
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix, prefix="xml2input_")
    tmp.write(raw)
    tmp.close()
    return Path(tmp.name), name, raw.decode("utf-8", errors="replace")


def _serialize_nuclide(n, name=None, text=None):
    betas = []
    for b in n["betas"]:
        betas.append({
            "from_e": b["from_e"],
            "to_e": b["to_e"],
            "q": b["q"],
            "intensity": b["intensity"],
            "charge": b.get("charge", -1),
            "z_target": b.get("z_target", 0),
            "a_target": b.get("a_target", 0),
        })
    gammas = []
    for src, lst in sorted(n["gammas"].items(), key=lambda kv: kv[0]):
        for g in lst:
            gammas.append({
                "from_e": g["from_e"],
                "to_e": g["to_e"],
                "q": g["q"],
                "intensity": g["intensity"],
            })
    return {
        "name": name or Path(n["path"]).name,
        "z": n["z"],
        "a": n["a"],
        "q_beta": n["q_beta"],
        "levels": list(n["levels"]),
        "betas": betas,
        "gammas": gammas,
        "text": text,
    }


def _inspect_path(path, display_name=None, text=None):
    if text is None:
        text = Path(path).read_text(encoding="utf-8", errors="replace")
    n = parse_nuclide_xml(path)
    return _serialize_nuclide(n, name=display_name, text=text)


def _paths_json(paths):
    out = []
    for i, p in enumerate(paths, start=1):
        out.append({
            "index": i,
            "weight": p.get("weight", 0.0),
            "feed": p["feed"],
            "chain": p["chain"],
            "beta_intensity": p["beta_intensity"],
            "gamma_intensities": p.get("gamma_intensities", []),
            "q_beta_branch": p.get("q_beta_branch"),
        })
    return out


def _translate_files(qbeta_path, gamma_path, names, e_max=None):
    result = run_translation([qbeta_path, gamma_path], e_max=e_max)
    arr = result["arr"]
    lim = float(max(abs(arr.min()), abs(arr.max())))
    if lim <= 0:
        lim = 1.0
    buf = io.BytesIO()
    title = f"{names[0]} + {names[1]}"
    plot_ladder(
        arr, buf, show=False,
        z_label="path probability",
        vmin=-lim, vmax=lim,
        title=title,
        verbose=False,
    )
    png_b64 = base64.b64encode(buf.getvalue()).decode("ascii")
    npy_buf = io.BytesIO()
    np.save(npy_buf, arr)
    npy_b64 = base64.b64encode(npy_buf.getvalue()).decode("ascii")
    npy_name = f"{Path(names[0]).stem}_{Path(names[1]).stem}_level_ladder_endpoint_normed.npy"
    return {
        "plot_png": png_b64,
        "npy": npy_b64,
        "npy_name": npy_name,
        "shape": list(arr.shape),
        "min": float(arr.min()),
        "max": float(arr.max()),
        "e_max": result["e_max"],
        "q_beta": result["q_beta"],
        "max_q_transition": result["max_q_transition"],
        "z_fermi": result["z_fermi"],
        "e_charge": result["e_charge"],
        "levels": result["canonical"],
        "paths": _paths_json(result["paths"]),
    }


@app.get("/")
def index():
    return send_from_directory(HERE, "index.html")


@app.get("/api/sample/<kind>")
def sample(kind):
    path = SAMPLE_QBETA if kind == "qbeta" else SAMPLE_GAMMA if kind == "gamma" else None
    if path is None or not path.is_file():
        return jsonify({"error": f"unknown sample {kind!r}"}), 404
    return jsonify(_inspect_path(path))


@app.post("/api/inspect")
def inspect():
    storage = request.files.get("file")
    if storage is None:
        return jsonify({"error": "missing file"}), 400
    tmp, name, text = _write_upload(storage)
    try:
        return jsonify(_inspect_path(tmp, display_name=name, text=text))
    except Exception as exc:
        return jsonify({"error": str(exc)}), 400
    finally:
        tmp.unlink(missing_ok=True)


@app.post("/api/translate")
def translate():
    e_max_raw = request.form.get("e_max", "").strip()
    e_max = float(e_max_raw) if e_max_raw else None
    use_sample = request.form.get("sample") == "1"

    tmp_paths = []
    try:
        if use_sample:
            q_path, g_path = SAMPLE_QBETA, SAMPLE_GAMMA
            names = [q_path.name, g_path.name]
        else:
            q_store = request.files.get("qbeta")
            g_store = request.files.get("gamma")
            if q_store is None or g_store is None:
                return jsonify({"error": "upload both a Qβ file and a γ file"}), 400
            q_path, q_name, _ = _write_upload(q_store)
            g_path, g_name, _ = _write_upload(g_store)
            tmp_paths.extend([q_path, g_path])
            names = [q_name, g_name]
        payload = _translate_files(q_path, g_path, names, e_max=e_max)
        return jsonify(payload)
    except Exception as exc:
        return jsonify({"error": str(exc)}), 400
    finally:
        for p in tmp_paths:
            p.unlink(missing_ok=True)


def main():
    app.run(host="127.0.0.1", port=5055, debug=False)


if __name__ == "__main__":
    main()
