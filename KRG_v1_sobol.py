# -*- coding: utf-8 -*-
# Kriging + Expected Improvement
# Pt = minimize
# Q = maximize
# 輸出：
# 1. EI suggested point（考慮不確定度）
# 2. Predicted optimal μ point（只看 mu，不考慮不確定度）
# 3. Best sample in data

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

from scipy.stats import norm, qmc
from scipy.optimize import differential_evolution


# =========================
# 路徑設定
# =========================
BASE_DIR = Path(__file__).parent

DATA_FILE = BASE_DIR / "LHS.txt"
MODEL_FILE = BASE_DIR / "krg_model.pkl"
RESULT_FILE = BASE_DIR / "ei_result.txt"


# =========================
# 設計變數上下界
# 依你的需求可自行修改
# 前 8 個變數會自動 round 成整數
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
# maximize=True: 最大化
# maximize=False: 最小化
# =========================
def expected_improvement(X_scaled, gp, f_best, maximize=True):
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
# 讀取資料
# 自動把欄位 rename 成 dh / Pt / Q
# =========================
def read_lhs_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"找不到檔案：{path}")

    df = pd.read_csv(path, sep=r"\s+|\t", engine="python")
    df.columns = df.columns.str.strip()

    rename_map = {}
    for c in df.columns:
        if "Δh" in c:
            rename_map[c] = "dh"
        if "ΔPt" in c:
            rename_map[c] = "Pt"
        if c.startswith("Q"):
            rename_map[c] = "Q"

    df = df.rename(columns=rename_map)

    required = [
        "r1", "r2", "r3", "r4",
        "m1", "m2", "m3", "m4",
        "dh", "Pt", "Q"
    ]

    for c in required:
        if c not in df.columns:
            raise RuntimeError(f"缺少欄位：{c}")

    for c in required:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna().reset_index(drop=True)
    return df


# =========================
# 訓練 Kriging
# =========================
def train_krg(X, y):
    scaler = MinMaxScaler()
    X_scaled = scaler.fit_transform(X)

    kernel = C(1.0, (1e-3, 1e5)) * RBF(
        length_scale=np.ones(X.shape[1]),
        length_scale_bounds=(1e-3, 1e5)
    )

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
# 把前 8 個變數變整數，並夾在 bounds 內
# =========================
def sanitize_design_vector(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float).copy()
    x = np.maximum(BOUNDS[:, 0], np.minimum(BOUNDS[:, 1], x))
    x[:8] = np.round(x[:8])
    x = np.maximum(BOUNDS[:, 0], np.minimum(BOUNDS[:, 1], x))
    return x


# =========================
# EI 搜尋點
# 用 Sobol 候選點找 EI 最大值
# =========================
def find_next_sample(gp, scaler, y_data, maximize):
    f_best = float(np.max(y_data) if maximize else np.min(y_data))

    sampler = qmc.Sobol(d=BOUNDS.shape[0], scramble=True, seed=0)
    U = sampler.random_base2(m=16)   # 65536 points

    xmin = BOUNDS[:, 0]
    xmax = BOUNDS[:, 1]

    X_rand = xmin + (xmax - xmin) * U
    X_rand[:, :8] = np.round(X_rand[:, :8])

    X_scaled = scaler.transform(X_rand)
    ei, mu, std = expected_improvement(X_scaled, gp, f_best, maximize=maximize)

    idx = int(np.argmax(ei))
    x_best = sanitize_design_vector(X_rand[idx])

    x_best_scaled = scaler.transform(x_best.reshape(1, -1))
    mu_best, std_best = gp.predict(x_best_scaled, return_std=True)

    return x_best, float(ei[idx]), float(mu_best[0]), float(std_best[0])


# =========================
# 真正的 surrogate μ 最佳點
# 不考慮不確定度，只看 μ
# Pt: argmin μ
# Q : argmax μ
# =========================
def find_true_mu_optimum(gp, scaler, maximize):
    bounds_list = [tuple(b) for b in BOUNDS]

    def objective(x):
        x = sanitize_design_vector(np.asarray(x, dtype=float))
        x_scaled = scaler.transform(x.reshape(1, -1))
        mu = float(gp.predict(x_scaled)[0])
        return -mu if maximize else mu

    result = differential_evolution(
        objective,
        bounds=bounds_list,
        strategy="best1bin",
        maxiter=200,
        popsize=15,
        tol=1e-6,
        mutation=(0.5, 1.0),
        recombination=0.7,
        polish=True,
        seed=0,
        updating="deferred",
        workers=1
    )

    x_opt = sanitize_design_vector(result.x)
    x_scaled = scaler.transform(x_opt.reshape(1, -1))
    mu_opt = float(gp.predict(x_scaled)[0])

    return x_opt, mu_opt


# =========================
# 主程式
# =========================
def main():
    df = read_lhs_table(DATA_FILE)

    X_cols = [
        "r1", "r2", "r3", "r4",
        "m1", "m2", "m3", "m4",
        "dh"
    ]

    print("選擇最佳化目標")
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

    # EI 點
    x_ei, ei_val, mu_ei, std_ei = find_next_sample(
        gp=gp,
        scaler=scaler,
        y_data=y,
        maximize=maximize
    )

    # 只看 μ 的最佳點
    x_mu, mu_best = find_true_mu_optimum(
        gp=gp,
        scaler=scaler,
        maximize=maximize
    )

    # 資料中最佳點
    best_idx = int(np.argmax(y) if maximize else np.argmin(y))

    lines = []
    lines.append("==== Kriging + Expected Improvement Result ====\n")
    lines.append(f"Optimization target = {y_col}")
    lines.append(f"Objective direction = {'maximize' if maximize else 'minimize'}")
    lines.append(f"Samples = {len(df)}")
    lines.append(f"Design variables = {X.shape[1]}")
    lines.append(f"Best data value = {y[best_idx]:.6f}\n")

    lines.append("Kernel:")
    lines.append(str(gp.kernel_) + "\n")

    lines.append("---- Bounds ----")
    for i, name in enumerate(X_cols):
        lines.append(f"{name}: [{BOUNDS[i,0]}, {BOUNDS[i,1]}]")
    lines.append("")

    lines.append("---- EI suggested point ----")
    for i, v in enumerate(X_cols):
        lines.append(f"{v} = {x_ei[i]:.6f}")
    lines.append(f"Predicted {y_col} = {mu_ei:.6f}")
    lines.append(f"Std = {std_ei:.6f}")
    lines.append(f"EI = {ei_val:.6f}\n")

    lines.append("---- Predicted optimal μ point ----")
    for i, v in enumerate(X_cols):
        lines.append(f"{v} = {x_mu[i]:.6f}")
    lines.append(f"Predicted best {y_col} = {mu_best:.6f}\n")

    lines.append("---- Best sample in data ----")
    if "Case" in df.columns:
        lines.append(f"Case = {int(df.loc[best_idx, 'Case'])}")
    for c in X_cols:
        lines.append(f"{c} = {df.loc[best_idx, c]}")
    lines.append(f"{y_col} = {df.loc[best_idx, y_col]:.6f}")

    RESULT_FILE.write_text("\n".join(lines), encoding="utf-8")

    bundle = {
        "gp": gp,
        "scaler": scaler,
        "X": X,
        "y": y,
        "target": y_col,
        "maximize": maximize,
        "bounds": BOUNDS,
        "x_ei": x_ei,
        "ei_val": ei_val,
        "mu_ei": mu_ei,
        "std_ei": std_ei,
        "x_mu": x_mu,
        "mu_best": mu_best,
    }
    joblib.dump(bundle, MODEL_FILE)


if __name__ == "__main__":
    main()
