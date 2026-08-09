# ruff: noqa
"""
apply_with_prior_single.py
==========================

Snakemake script: PyTorch inference for ONE AFID.

Called by the ``applyfidmodel_single`` rule which has wildcard ``{afid}``
(zero-padded two digits: 01–32).  Snakemake runs ``--cores N`` instances of
this rule concurrently — that is the parallelism mechanism.  No Python
multiprocessing is used here.

Snakemake I/O
-------------
    input:
        t1w   – preprocessed NIfTI (.nii.gz)
        prior – per-subject MNI-registered FCSV with 32 prior locations
    output:
        coord – plain-text file: "x y z\\n" for this AFID
    wildcards:
        afid  – zero-padded AFID number, e.g. "01"
"""

import os

# Set environment variables BEFORE any other imports to ensure they take effect
if "snakemake" in globals() and hasattr(snakemake, "threads"):
    n_threads = str(snakemake.threads)
    os.environ["OMP_NUM_THREADS"] = n_threads
    os.environ["MKL_NUM_THREADS"] = n_threads
    os.environ["OPENBLAS_NUM_THREADS"] = n_threads
    os.environ["VECLIB_MAXIMUM_THREADS"] = n_threads
    os.environ["NUMEXPR_NUM_THREADS"] = n_threads

import time
import warnings
from pathlib import Path
from typing import List, Optional, Tuple

import nibabel as nib
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from numpy.typing import NDArray

if "snakemake" in globals() and hasattr(snakemake, "threads"):
    # print(f"DEBUG: Snakemake threads = {snakemake.threads}")
    # print(f"DEBUG: torch.get_num_threads() (before) = {torch.get_num_threads()}")
    torch.set_num_threads(snakemake.threads)
    # print(f"DEBUG: torch.get_num_threads() (after) = {torch.get_num_threads()}")

warnings.filterwarnings("ignore")


# ===================================================================
# FCSV helpers
# ===================================================================


def load_fcsv(fcsv_path) -> pd.DataFrame:
    return pd.read_csv(fcsv_path, sep=",", header=2)


def get_fid(fcsv_df: pd.DataFrame, fid_label: int) -> NDArray:
    return fcsv_df.loc[fid_label - 1, ["x", "y", "z"]].to_numpy(
        dtype="single", copy=True
    )


def fid_world2voxel(fid_world: NDArray, nii_affine: NDArray) -> NDArray:
    inv = np.linalg.inv(nii_affine)
    return np.rint(inv[:3, :3].dot(fid_world) + inv[:3, 3]).astype(int)


def fid_voxel2world(fid_voxel: NDArray, nii_affine: NDArray) -> NDArray:
    return (nii_affine[:3, :3].dot(fid_voxel) + nii_affine[:3, 3]).astype(float)


# ===================================================================
# nnUNet architecture — names MUST match checkpoint keys exactly
# (encoder_blocks, decoder_blocks, deep_supervision_heads)
# ===================================================================


class ConvBlock(nn.Module):
    def __init__(self, in_channels: int, out_channels: int) -> None:
        super().__init__()
        self.conv_block = nn.Sequential(
            nn.Conv3d(in_channels, out_channels, kernel_size=3, padding=1, bias=True),
            nn.InstanceNorm3d(out_channels, affine=True),
            nn.LeakyReLU(negative_slope=0.01, inplace=True),
            nn.Conv3d(out_channels, out_channels, kernel_size=3, padding=1, bias=True),
            nn.InstanceNorm3d(out_channels, affine=True),
            nn.LeakyReLU(negative_slope=0.01, inplace=True),
        )

    def forward(self, x):
        return self.conv_block(x)


class EncoderBlock(nn.Module):
    def __init__(self, in_channels: int, out_channels: int) -> None:
        super().__init__()
        self.conv_block = ConvBlock(in_channels, out_channels)
        self.downsample = nn.Conv3d(
            out_channels, out_channels, kernel_size=3, stride=2, padding=1
        )

    def forward(self, x):
        skip = self.conv_block(x)
        return skip, self.downsample(skip)


class DecoderBlock(nn.Module):
    def __init__(self, in_channels: int, out_channels: int) -> None:
        super().__init__()
        self.upsample = nn.ConvTranspose3d(
            in_channels, in_channels // 2, kernel_size=2, stride=2
        )
        self.conv_block = ConvBlock(in_channels, out_channels)

    def forward(self, x, skip):
        return self.conv_block(torch.cat([self.upsample(x), skip], 1))


class nnUNet(nn.Module):
    def __init__(
        self, in_channels: int, out_channels: int, features: Optional[List[int]] = None
    ) -> None:
        super().__init__()
        if features is None:
            features = [32, 64, 128, 256, 320]
        self.encoder_blocks = nn.ModuleList()
        prev_ch = in_channels
        for feat in features[:-1]:
            self.encoder_blocks.append(EncoderBlock(prev_ch, feat))
            prev_ch = feat
        self.bottleneck = ConvBlock(features[-2], features[-1])
        self.decoder_blocks = nn.ModuleList()
        rev = list(reversed(features))
        for i in range(len(rev) - 1):
            self.decoder_blocks.append(DecoderBlock(rev[i], rev[i + 1]))
        self.deep_supervision_heads = nn.ModuleList()
        for feat in reversed(features[:-1]):
            self.deep_supervision_heads.append(
                nn.Conv3d(feat, out_channels, kernel_size=1)
            )

    def forward(self, x):
        skips = []
        for enc in self.encoder_blocks:
            skip, x = enc(x)
            skips.append(skip)
        x = self.bottleneck(x)
        skips = list(reversed(skips))
        for i, dec in enumerate(self.decoder_blocks):
            x = dec(x, skips[i])
        return self.deep_supervision_heads[-1](x)


class nnUNet_VanillaUNet(nnUNet):
    """Thin wrapper matching the Lightning checkpoint's class."""

    def __init__(
        self,
        in_channels: int = 1,
        out_channels: int = 1,
        features: Optional[List[int]] = None,
    ) -> None:
        super().__init__(in_channels, out_channels, features)


def _load_model(
    ckpt_path: str, features: List[int], device: torch.device
) -> nnUNet_VanillaUNet:
    model = nnUNet_VanillaUNet(in_channels=1, out_channels=1, features=features)
    ckpt = torch.load(ckpt_path, map_location="cpu")
    raw_sd = ckpt.get("state_dict", ckpt)
    cleaned_sd = {
        (k.replace("model.", "", 1) if k.startswith("model.") else k): v
        for k, v in raw_sd.items()
    }
    model.load_state_dict(cleaned_sd, strict=False)
    return model.to(device).eval()


# ===================================================================
# Single-AFID inference
# ===================================================================


def infer_single_afid(
    fid: int,
    ckpt_path: str,
    img: nib.nifti1.Nifti1Image,
    prior_fcsv_path,
    patch_size: int = 64,
    batch_size: int = 5,
    device_str: str = "cpu",
    features: Optional[List[int]] = None,
) -> Tuple[NDArray, torch.Tensor]:
    """Run 5-patch inference for a single AFID.

    Returns
    -------
    pred_world : NDArray  – predicted world coords [x, y, z]
    prob_map   : Tensor  – full-volume Gaussian-weighted probability map [D,H,W]
    """
    if features is None:
        features = [16, 32, 64]

    device = torch.device(
        device_str if (device_str == "cpu" or torch.cuda.is_available()) else "cpu"
    )

    # ---- Image ----
    image_np = img.get_fdata().astype(np.float32)
    if not np.isfinite(image_np).all():
        finite = np.isfinite(image_np)
        lo = image_np[finite].min() if finite.any() else 0.0
        hi = image_np[finite].max() if finite.any() else 1.0
        image_np = np.nan_to_num(image_np, nan=0.0, posinf=hi, neginf=lo)
    image_tensor = torch.from_numpy(image_np)
    affine = img.affine
    vol_shape = image_tensor.shape
    D, H, W = vol_shape

    # ---- Prior centre ----
    prior_df = load_fcsv(prior_fcsv_path)
    prior_world = get_fid(prior_df, fid)
    prior_vox = fid_world2voxel(prior_world, affine)
    cz = int(np.clip(prior_vox[0], 0, D - 1))
    cy = int(np.clip(prior_vox[1], 0, H - 1))
    cx = int(np.clip(prior_vox[2], 0, W - 1))

    # ---- Load model ----
    t0 = time.perf_counter()
    model = _load_model(ckpt_path, features, device)
    t_load = time.perf_counter() - t0

    # ---- Gaussian map ----
    ps = patch_size

    def _g1d(n):
        s = n * 0.125
        c = n // 2
        x = torch.arange(n, dtype=torch.float32)
        return torch.exp(-((x - c) ** 2) / (2 * s**2))

    gmap = torch.einsum("z,y,x->zyx", _g1d(ps), _g1d(ps), _g1d(ps))
    gmap = gmap / gmap.max()

    # ---- 7 patches: centre + ±x + ±y + ±z ----
    half = ps // 2
    offsets = [
        (0, 0, 0),  # centre
        (0, 0, -half),  # left   (x−)
        (0, 0, +half),  # right  (x+)
        (0, -half, 0),  # up     (y−)
        (0, +half, 0),  # down   (y+)
        (-half, 0, 0),  # in     (z−)
        (+half, 0, 0),  # out    (z+)
    ]
    patch_coords = []
    for dz, dy, dx in offsets:
        zs = max(0, min(cz + dz - half, D - ps))
        ys = max(0, min(cy + dy - half, H - ps))
        xs = max(0, min(cx + dx - half, W - ps))
        patch_coords.append(
            (slice(zs, zs + ps), slice(ys, ys + ps), slice(xs, xs + ps))
        )

    # ---- Extract & normalise ----
    t0 = time.perf_counter()
    model_inputs = []
    for slcs in patch_coords:
        patch = image_tensor[slcs[0], slcs[1], slcs[2]]
        pz = max(0, ps - patch.shape[0])
        py = max(0, ps - patch.shape[1])
        px_p = max(0, ps - patch.shape[2])
        if pz or py or px_p:
            patch = F.pad(patch, (0, px_p, 0, py, 0, pz))
        lo, hi = patch.min(), patch.max()
        patch = (patch - lo) / (hi - lo + 1e-8) if hi > lo else patch
        model_inputs.append(patch.unsqueeze(0))
    t_patch = time.perf_counter() - t0

    # ---- Inference ----
    t0 = time.perf_counter()
    cpu0 = time.process_time()
    predictions = []
    for i in range(0, len(model_inputs), batch_size):
        chunk = torch.stack(model_inputs[i : i + batch_size]).to(device)
        with torch.inference_mode():
            preds = model(chunk)
        predictions.extend(p.cpu() for p in preds)
        del chunk, preds
    t_infer = time.perf_counter() - t0
    cpu_infer = time.process_time() - cpu0

    # ---- Reconstruct ----
    t0 = time.perf_counter()
    out = torch.zeros(vol_shape, dtype=torch.float32)
    cnt = torch.zeros(vol_shape, dtype=torch.float32)
    for pred, (zs, ys, xs) in zip(predictions, patch_coords):
        out[zs, ys, xs] += pred[0] * gmap
        cnt[zs, ys, xs] += gmap
    out /= cnt.clamp(min=1e-8)
    t_recon = time.perf_counter() - t0

    # ---- Peak → world ----
    idx = np.unravel_index(np.argmax(out.numpy()), vol_shape)
    pred_vox = tuple(int(c) for c in idx)
    pred_world = fid_voxel2world(np.array(pred_vox, dtype=float), affine)

    t_total = t_load + t_patch + t_infer + t_recon
    print(
        f"  [OK]  AFID {fid:02d}  "
        f"load={t_load:.2f}s  patch={t_patch:.2f}s  "
        f"infer={t_infer:.2f}s (cpu={cpu_infer:.2f}s)  "
        f"recon={t_recon:.2f}s  total={t_total:.2f}s  |  "
        f"pred_vox={pred_vox}  "
        f"world=[{pred_world[0]:.1f}, {pred_world[1]:.1f}, {pred_world[2]:.1f}]  "
        f"[dev={device}, threads={torch.get_num_threads()}]"
    )
    return pred_world, out  # out = full-volume probability map [D,H,W]


def _infer_single_from_state(
    fid: int,
    image_tensor: torch.Tensor,
    affine: NDArray,
    prior_df: pd.DataFrame,
    model: nn.Module,
    gmap: torch.Tensor,
    patch_size: int,
    batch_size: int,
    device: torch.device,
) -> Tuple[NDArray, torch.Tensor]:
    """Core inference that reuses precomputed image tensor, prior, model, and gmap."""
    vol_shape = tuple(int(d) for d in image_tensor.shape)
    D, H, W = vol_shape

    # ---- Prior centre ----
    prior_world = get_fid(prior_df, fid)
    prior_vox = fid_world2voxel(prior_world, affine)
    cz = int(np.clip(prior_vox[0], 0, D - 1))
    cy = int(np.clip(prior_vox[1], 0, H - 1))
    cx = int(np.clip(prior_vox[2], 0, W - 1))

    # ---- 7 patches: centre + ±x + ±y + ±z ----
    ps = patch_size
    half = ps // 2
    offsets = [
        (0, 0, 0),
        (0, 0, -half),
        (0, 0, +half),
        (0, -half, 0),
        (0, +half, 0),
        (-half, 0, 0),
        (+half, 0, 0),
    ]
    patch_coords = []
    for dz, dy, dx in offsets:
        zs = max(0, min(cz + dz - half, D - ps))
        ys = max(0, min(cy + dy - half, H - ps))
        xs = max(0, min(cx + dx - half, W - ps))
        patch_coords.append(
            (slice(zs, zs + ps), slice(ys, ys + ps), slice(xs, xs + ps))
        )

    # ---- Extract & normalise ----
    t0 = time.perf_counter()
    model_inputs = []
    for slcs in patch_coords:
        patch = image_tensor[slcs[0], slcs[1], slcs[2]]
        pz = max(0, ps - patch.shape[0])
        py = max(0, ps - patch.shape[1])
        px_p = max(0, ps - patch.shape[2])
        if pz or py or px_p:
            patch = F.pad(patch, (0, px_p, 0, py, 0, pz))
        lo, hi = patch.min(), patch.max()
        patch = (patch - lo) / (hi - lo + 1e-8) if hi > lo else patch
        model_inputs.append(patch.unsqueeze(0))
    t_patch = time.perf_counter() - t0

    # ---- Inference ----
    t0 = time.perf_counter()
    cpu0 = time.process_time()
    predictions = []
    for i in range(0, len(model_inputs), batch_size):
        chunk = torch.stack(model_inputs[i : i + batch_size]).to(device)
        with torch.inference_mode():
            preds = model(chunk)
        predictions.extend(p.cpu() for p in preds)
        del chunk, preds
    t_infer = time.perf_counter() - t0
    cpu_infer = time.process_time() - cpu0

    # ---- Reconstruct ----
    t0 = time.perf_counter()
    out = torch.zeros(vol_shape, dtype=torch.float32)
    cnt = torch.zeros(vol_shape, dtype=torch.float32)
    for pred, (zs, ys, xs) in zip(predictions, patch_coords):
        out[zs, ys, xs] += pred[0] * gmap
        cnt[zs, ys, xs] += gmap
    out /= cnt.clamp(min=1e-8)
    t_recon = time.perf_counter() - t0

    # ---- Peak → world ----
    idx = np.unravel_index(np.argmax(out.numpy()), vol_shape)
    pred_vox = tuple(int(c) for c in idx)
    pred_world = fid_voxel2world(np.array(pred_vox, dtype=float), affine)

    t_total = t_patch + t_infer + t_recon
    print(
        f"  [OK]  AFID {fid:02d}  "
        f"patch={t_patch:.2f}s  "
        f"infer={t_infer:.2f}s (cpu={cpu_infer:.2f}s)  "
        f"recon={t_recon:.2f}s  total={t_total:.2f}s  |  "
        f"pred_vox={pred_vox}  "
        f"world=[{pred_world[0]:.1f}, {pred_world[1]:.1f}, {pred_world[2]:.1f}]  "
        f"[dev={device}, threads={torch.get_num_threads()}]"
    )
    return pred_world, out


# ===================================================================
# Snakemake entry point (All 32 AFIDs)
# ===================================================================

import os
from pathlib import Path
import csv

AFIDS_FIELDNAMES = [
    "id",
    "x",
    "y",
    "z",
    "ow",
    "ox",
    "oy",
    "oz",
    "vis",
    "sel",
    "lock",
    "label",
    "desc",
    "associatedNodeID",
]

AFID_DESCRIPTIONS = [
    "AC",
    "PC",
    "Infracollicular Sulcus",
    "PMJ",
    "Superior IPF",
    "Right Superior LMS",
    "Left Superior LMS",
    "Right Inferior LMS",
    "Left Inferior LMS",
    "Culmen",
    "Intermammillary Sulcus",
    "Right Mammilary Body",
    "Left Mammilary Body",
    "Pineal Gland",
    "Right LV at AC",
    "Left LV at AC",
    "Right LV at PC",
    "Left LV at PC",
    "Genu of CC",
    "Splenium of CC",
    "Right AL Temporal Horn",
    "Left AL Tempral Horn",
    "R. Sup. AM Temporal Horn",
    "L. Sup. AM Temporal Horn",
    "R Inf. AM Temp Horn",
    "L Inf. AM Temp Horn",
    "Right IG Origin",
    "Left IG Origin",
    "R Ventral Occipital Horn",
    "L Ventral Occipital Horn",
    "R Olfactory Fundus",
    "L Olfactory Fundus",
]


def write_combined_fcsv(afid_coords, fcsv_output):
    header_lines = [
        "# Markups fiducial file version = 4.10\n",
        "# CoordinateSystem = 0\n",
        "# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n",
    ]
    rows = []
    for lbl in range(1, 33):
        c = afid_coords.get(lbl, np.array([0.0, 0.0, 0.0]))
        desc = AFID_DESCRIPTIONS[lbl - 1] if lbl <= len(AFID_DESCRIPTIONS) else ""
        rows.append(
            {
                "id": lbl,
                "x": c[0],
                "y": c[1],
                "z": c[2],
                "ow": "0.000",
                "ox": "0.000",
                "oy": "0.000",
                "oz": "1.000",
                "vis": 1,
                "sel": 1,
                "lock": 1,
                "label": lbl,
                "desc": desc,
                "associatedNodeID": "",
            }
        )
    with Path(fcsv_output).open("w", encoding="utf-8", newline="") as f:
        for line in header_lines:
            f.write(line)
        writer = csv.DictWriter(f, fieldnames=AFIDS_FIELDNAMES)
        for row in rows:
            writer.writerow(row)


cfg_block = snakemake.config.get("afids_inference", {})
if not cfg_block:
    raise ValueError("Missing 'afids_inference' block in snakebids.yml.")

checkpoints = snakemake.params.ckpts
img = nib.nifti1.load(snakemake.input.t1w)
prior = snakemake.input.prior

# Precompute image tensor and clean NaNs/Infs once per subject
image_np = img.get_fdata().astype(np.float32)
if not np.isfinite(image_np).all():
    finite = np.isfinite(image_np)
    lo = image_np[finite].min() if finite.any() else 0.0
    hi = image_np[finite].max() if finite.any() else 1.0
    image_np = np.nan_to_num(image_np, nan=0.0, posinf=hi, neginf=lo)
image_tensor = torch.from_numpy(image_np)
affine = img.affine

# Load prior FCSV once
prior_df = load_fcsv(prior)

# Precompute Gaussian map once
ps = int(cfg_block.get("patch_size", 64))


def _g1d_all(n: int) -> torch.Tensor:
    s = n * 0.125
    c = n // 2
    x = torch.arange(n, dtype=torch.float32)
    return torch.exp(-((x - c) ** 2) / (2 * s**2))


gmap = torch.einsum("z,y,x->zyx", _g1d_all(ps), _g1d_all(ps), _g1d_all(ps))
gmap = gmap / gmap.max()

batch_size = snakemake.config.get("inference_batch_size") or 7
device_str = cfg_block.get("device", "cpu")
device = torch.device(
    device_str if (device_str == "cpu" or torch.cuda.is_available()) else "cpu"
)

# Optional: cache models by checkpoint path
features = cfg_block.get("features", [16, 32, 64])
model_cache = {}

combined_coords = {}

for fid in range(1, 33):
    key = f"afid_{fid:02d}"
    ckpt_path = checkpoints.get(key)
    out_coord = snakemake.output.coords[fid - 1]
    out_prob = snakemake.output.probs[fid - 1]

    if not ckpt_path or not Path(ckpt_path).exists():
        print(f"  [SKIP] AFID {fid:02d} checkpoint not found: {ckpt_path}")
        combined_coords[fid] = np.array([0.0, 0.0, 0.0])
        Path(out_coord).parent.mkdir(parents=True, exist_ok=True)
        with open(out_coord, "w") as f:
            f.write("0.0 0.0 0.0\n")
        continue

    # Load or reuse model for this checkpoint
    model = model_cache.get(ckpt_path)
    if model is None:
        model = _load_model(ckpt_path, features, device)
        model_cache[ckpt_path] = model

    pred_world, prob_map = _infer_single_from_state(
        fid=fid,
        image_tensor=image_tensor,
        affine=affine,
        prior_df=prior_df,
        model=model,
        gmap=gmap,
        patch_size=ps,
        batch_size=batch_size,
        device=device,
    )
    combined_coords[fid] = pred_world

    Path(out_coord).parent.mkdir(parents=True, exist_ok=True)
    with open(out_coord, "w") as f:
        f.write(f"{pred_world[0]} {pred_world[1]} {pred_world[2]}\n")

    Path(out_prob).parent.mkdir(parents=True, exist_ok=True)
    nib.save(
        nib.Nifti1Image(prob_map.numpy(), affine, img.header),
        out_prob,
    )

write_combined_fcsv(combined_coords, snakemake.output.fcsv)
