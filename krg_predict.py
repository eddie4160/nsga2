import joblib
import numpy as np
from pathlib import Path


_CACHE = {}


def _load(model_filename):
    model_path = Path(__file__).parent / model_filename
    if model_filename not in _CACHE:
        _CACHE[model_filename] = joblib.load(model_path)
    return _CACHE[model_filename]


def _get_target_scaler(pack):
    for key in ("y_scaler", "target_scaler", "output_scaler"):
        if key in pack:
            return pack[key]
    return None


def _inverse_transform_mu_and_std(target_scaler, mu, std):
    mu = np.asarray(mu, dtype=float).reshape(-1, 1)
    std = np.asarray(std, dtype=float).reshape(-1)
    if target_scaler is None:
        return mu.ravel(), std
    mu = target_scaler.inverse_transform(mu).ravel()
    scale = getattr(target_scaler, "scale_", None)
    if scale is not None:
        std = std * float(np.ravel(scale)[0])
    return mu, std


def predict(model_filename, values):
    pack = _load(model_filename)
    x = np.asarray(values, dtype=float).reshape(1, -1)
    scaler = pack.get("scaler")
    gp = pack.get("gp", pack.get("model"))
    if gp is None:
        raise KeyError("Expected 'gp' or 'model' in surrogate package.")
    if scaler is not None:
        x = scaler.transform(x)
    mu, std = gp.predict(x, return_std=True)
    mu, std = _inverse_transform_mu_and_std(_get_target_scaler(pack), mu, std)
    return float(mu[0]), float(std[0])
