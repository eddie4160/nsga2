# -*- coding: utf-8 -*-
# Kriging + EI optimizer
# integer grid + dh coarse-to-fine search
# 輸出：
# 1) EI suggested point
# 2) Predicted optimal μ point
# 3) Best sample in data

from __future__ import annotations

import numpy as np
import pandas as pd
from pathlib import Path
import joblib
import warnings

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
from sklearn.preprocessing import MinMaxScaler
from sklearn.exceptions import ConvergenceWarning

from scipy.stats import norm


BASE_DIR = Path(__file__).parent

DATA_FILE = BASE_DIR / "LHS.txt"
MODEL_FILE = BASE_DIR / "krg_model.pkl"
RESULT_FILE = BASE_DIR / "ei_result.txt"


# =========================
# bounds
# =========================
BOUNDS = np.array([
    [1, 20],   # r1
    [1, 20],   # r2
    [1, 20],   # r3
    [1, 20],   # r4
    [1, 6],    # m1
    [1, 6],    # m2
    [1, 6],    # m3
    [1, 6],    # m4
    [-5, 20],  # dh
], dtype=float)


# =========================
# EI
# maximize=True  -> for Q
# maximize=False -> for Pt
# =========================
def expected_improvement(X_scaled, gp, f_best, maximize):
    mu, std = gp.predict(X_scaled, return_std=True)
    std = np.maximum(std, 1e-12)

    if maximize:
        imp = mu - f_best
    else:
        imp = f_best - mu

    z = imp / std
    ei = imp * norm.cdf(z) + std * norm.pdf(z)
    ei[std < 1e-12] = 0.0

    return ei, mu, std


# =========================
# read data
# =========================
def read_lhs_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"找不到檔案：{path}")

    df = pd.read_csv(path, sep=r"\s+|\t", engine="python")
    df.columns = df.columns.str.strip()

    rename = {}
    for c in df.columns:
        if "Δh" in c:
            rename[c] = "dh"
        if "ΔPt" in c:
            rename[c] = "Pt"
        if c.startswith("Q"):
            rename[c] = "Q"

    df = df.rename(columns=rename)

    required_cols = [
        "r1", "r2", "r3", "r4",
        "m1", "m2", "m3", "m4",
        "dh", "Pt", "Q"
    ]

    for c in required_cols:
        if c not in df.columns:
            raise RuntimeError(f"缺少欄位: {c}")
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna().reset_index(drop=True)
    return df


# =========================
# train kriging
# =========================
def train_krg(X, y):
    scaler = MinMaxScaler()
    X_scaled = scaler.fit_transform(X)

    kernel = C(1.0, (1e-3, 1e5)) * RBF(length_scale=np.ones(X.shape[1]))

    gp = GaussianProcessRegressor(
        kernel=kernel,
        alpha=1e-6,
        normalize_y=True,
        n_restarts_optimizer=20,
        random_state=0
    )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        gp.fit(X_scaled, y)

    return gp, scaler


# =========================
# integer grid samples
# 隨機抽整數組合
# =========================
def integer_samples(n=20000):
    rng = np.random.default_rng(0)

    r = rng.integers(1, 21, (n, 4))
    m = rng.integers(1, 7,  (n, 4))

    return np.hstack([r, m]).astype(float)


# =========================
# dh coarse-to-fine refine
# 對固定整數組合 x_int
# 同時找：
# 1) local EI best
# 2) local mu best
# =========================
def refine_dh(gp, scaler, x_int, maximize, y_best):
    dh_min, dh_max = BOUNDS[8]

    local_ei_best = None
    local_mu_best = None

    # 四階段 coarse-to-fine
    stages = [
        (0.1,    None),
        (0.01,   0.2),
        (0.001,  0.02),
        (0.00001, 0.002),
    ]

    # EI 用自己的中心，mu 用自己的中心
    best_dh_ei = None
    best_ei = None
    best_mu_at_ei = None
    best_std_at_ei = None

    best_dh_mu = None
    best_mu = None
    best_std_at_mu = None
    best_ei_at_mu = None

    for step, window in stages:
        # ---- EI branch ----
        if window is None or best_dh_ei is None:
            dh_list_ei = np.arange(dh_min, dh_max + step, step)
        else:
            lo = max(dh_min, best_dh_ei - window)
            hi = min(dh_max, best_dh_ei + window)
            dh_list_ei = np.arange(lo, hi + step, step)

        X_ei = np.column_stack([
            np.repeat(x_int[0], len(dh_list_ei)),
            np.repeat(x_int[1], len(dh_list_ei)),
            np.repeat(x_int[2], len(dh_list_ei)),
            np.repeat(x_int[3], len(dh_list_ei)),
            np.repeat(x_int[4], len(dh_list_ei)),
            np.repeat(x_int[5], len(dh_list_ei)),
            np.repeat(x_int[6], len(dh_list_ei)),
            np.repeat(x_int[7], len(dh_list_ei)),
            dh_list_ei
        ])

        X_scaled_ei = scaler.transform(X_ei)
        ei_vals_ei, mu_vals_ei, std_vals_ei = expected_improvement(
            X_scaled_ei, gp, y_best, maximize
        )

        idx_ei = int(np.argmax(ei_vals_ei))
        best_dh_ei = float(dh_list_ei[idx_ei])
        best_ei = float(ei_vals_ei[idx_ei])
        best_mu_at_ei = float(mu_vals_ei[idx_ei])
        best_std_at_ei = float(std_vals_ei[idx_ei])

        # ---- mu branch ----
        if window is None or best_dh_mu is None:
            dh_list_mu = np.arange(dh_min, dh_max + step, step)
        else:
            lo = max(dh_min, best_dh_mu - window)
            hi = min(dh_max, best_dh_mu + window)
            dh_list_mu = np.arange(lo, hi + step, step)

        X_mu = np.column_stack([
            np.repeat(x_int[0], len(dh_list_mu)),
            np.repeat(x_int[1], len(dh_list_mu)),
            np.repeat(x_int[2], len(dh_list_mu)),
            np.repeat(x_int[3], len(dh_list_mu)),
            np.repeat(x_int[4], len(dh_list_mu)),
            np.repeat(x_int[5], len(dh_list_mu)),
            np.repeat(x_int[6], len(dh_list_mu)),
            np.repeat(x_int[7], len(dh_list_mu)),
            dh_list_mu
        ])

        X_scaled_mu = scaler.transform(X_mu)
        ei_vals_mu, mu_vals_mu, std_vals_mu = expected_improvement(
            X_scaled_mu, gp, y_best, maximize
        )

        idx_mu = int(np.argmax(mu_vals_mu) if maximize else np.argmin(mu_vals_mu))
        best_dh_mu = float(dh_list_mu[idx_mu])
        best_mu = float(mu_vals_mu[idx_mu])
        best_std_at_mu = float(std_vals_mu[idx_mu])
        best_ei_at_mu = float(ei_vals_mu[idx_mu])

    local_ei_best = {
        "x": np.append(x_int, best_dh_ei),
        "mu": best_mu_at_ei,
        "std": best_std_at_ei,
        "ei": best_ei
    }

    local_mu_best = {
        "x": np.append(x_int, best_dh_mu),
        "mu": best_mu,
        "std": best_std_at_mu,
        "ei": best_ei_at_mu
    }

    return local_ei_best, local_mu_best


# =========================
# main search
# 全域找：
# 1) EI suggested point
# 2) Predicted optimal μ point
# =========================
def search_space(gp, scaler, y, maximize):
    y_best = float(np.max(y) if maximize else np.min(y))
    X_int = integer_samples()

    global_ei_best = None
    global_mu_best = None

    for x_int in X_int:
        local_ei_best, local_mu_best = refine_dh(
            gp, scaler, x_int, maximize, y_best
        )

        # ---- update global EI best ----
        if global_ei_best is None or local_ei_best["ei"] > global_ei_best["ei"]:
            global_ei_best = local_ei_best

        # ---- update global mu best ----
        if global_mu_best is None:
            global_mu_best = local_mu_best
        else:
            if maximize:
                if local_mu_best["mu"] > global_mu_best["mu"]:
                    global_mu_best = local_mu_best
            else:
                if local_mu_best["mu"] < global_mu_best["mu"]:
                    global_mu_best = local_mu_best

    return global_ei_best, global_mu_best


# =========================
# main
# =========================
def main():
    df = read_lhs_table(DATA_FILE)

    X_cols = [
        "r1", "r2", "r3", "r4",
        "m1", "m2", "m3", "m4",
        "dh"
    ]

    print("1 = ΔPt (min)")
    print("2 = Q (max)")
    choice = input("請輸入 1 或 2 : ").strip()

    if choice == "2":
        y_col = "Q"
        maximize = True
    else:
        y_col = "Pt"
        maximize = False

    X = df[X_cols].to_numpy(float)
    y = df[y_col].to_numpy(float)

    gp, scaler = train_krg(X, y)

    ei_best, mu_best = search_space(gp, scaler, y, maximize)

    best_idx = int(np.argmax(y) if maximize else np.argmin(y))

    lines = []
    lines.append("==== Kriging + Expected Improvement Result ====\n")
    lines.append(f"Optimization target = {y_col}")
    lines.append(f"Samples = {len(df)}")
    lines.append(f"Best data value = {y[best_idx]:.6f}\n")

    lines.append("---- EI suggested point ----")
    for i, v in enumerate(X_cols):
        lines.append(f"{v} = {ei_best['x'][i]:.6f}")
    lines.append(f"Predicted {y_col} = {ei_best['mu']:.6f}")
    lines.append(f"Std = {ei_best['std']:.6f}")
    lines.append(f"EI = {ei_best['ei']:.6f}\n")

    lines.append("---- Predicted optimal μ point ----")
    for i, v in enumerate(X_cols):
        lines.append(f"{v} = {mu_best['x'][i]:.6f}")
    lines.append(f"Predicted best {y_col} = {mu_best['mu']:.6f}")
    lines.append(f"Std = {mu_best['std']:.6f}")
    lines.append(f"EI = {mu_best['ei']:.6f}\n")

    lines.append("---- Best sample in data ----")
    if "Case" in df.columns:
        lines.append(f"Case = {int(df.loc[best_idx, 'Case'])}")

    for c in X_cols:
        lines.append(f"{c} = {df.loc[best_idx, c]}")

    lines.append(f"{y_col} = {df.loc[best_idx, y_col]:.6f}")

    RESULT_FILE.write_text("\n".join(lines), encoding="utf-8")

    joblib.dump(
        {
            "gp": gp,
            "scaler": scaler,
            "target": y_col,
            "maximize": maximize,
            "bounds": BOUNDS,
            "ei_best": ei_best,
            "mu_best": mu_best,
        },
        MODEL_FILE
    )


if __name__ == "__main__":
    main()