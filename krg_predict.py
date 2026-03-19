import joblib
import numpy as np
from pathlib import Path


_CACHE = {}


def _load(model_filename):
    model_path = Path(__file__).parent / model_filename
    if model_filename not in _CACHE:
        _CACHE[model_filename] = joblib.load(model_path)
    return _CACHE[model_filename]


def predict(model_filename, values):
    pack = _load(model_filename)
    x = np.asarray(values, dtype=float).reshape(1, -1)
    x_scaled = pack["scaler"].transform(x)
    mu, std = pack["gp"].predict(x_scaled, return_std=True)
    return float(mu[0]), float(std[0])
