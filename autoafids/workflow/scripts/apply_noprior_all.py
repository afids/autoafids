# ruff: noqa
"""
apply_noprior_single.py
=======================

Snakemake script: whole-volume sliding-window PyTorch inference for ONE AFID.
Does NOT require a prior FCSV — the model scans the entire image.

Called by the ``applyfidmodel_noprior_single`` rule which has wildcard ``{afid}``.
Snakemake runs ``--cores N`` instances concurrently — that is the parallelism.

Snakemake I/O
-------------
    input:
        t1w   – preprocessed NIfTI (.nii.gz)   [no prior needed]
    output:
        coord – plain-text file: "x y z\\n" for this AFID
        prob  – full-volume probability map as NIfTI (.nii.gz)
    wildcards:
        afid  – zero-padded AFID number, e.g. "01"
"""

import time
import warnings
from pathlib import Path
from typing import List, Optional, Tuple

import nibabel as nib
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

import os

# Each Snakemake job = 1 process. Limit PyTorch threads so that
# --cores N means N truly parallel jobs without over-subscription thrashing.
if "snakemake" in globals() and hasattr(snakemake, "threads"):
    n_threads = str(snakemake.threads)
    os.environ["OMP_NUM_THREADS"] = n_threads
    os.environ["MKL_NUM_THREADS"] = n_threads
    os.environ["OPENBLAS_NUM_THREADS"] = n_threads
    os.environ["VECLIB_MAXIMUM_THREADS"] = n_threads
    os.environ["NUMEXPR_NUM_THREADS"] = n_threads
    torch.set_num_threads(snakemake.threads)

from numpy.typing import NDArray

warnings.filterwarnings("ignore")


# ===================================================================
# FCSV helpers
# ===================================================================

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
    def forward(self, x): return self.conv_block(x)


class EncoderBlock(nn.Module):
    def __init__(self, in_channels: int, out_channels: int) -> None:
        super().__init__()
        self.conv_block = ConvBlock(in_channels, out_channels)
        self.downsample = nn.Conv3d(out_channels, out_channels, kernel_size=3, stride=2, padding=1)
    def forward(self, x):
        skip = self.conv_block(x)
        return skip, self.downsample(skip)


class DecoderBlock(nn.Module):
    def __init__(self, in_channels: int, out_channels: int) -> None:
        super().__init__()
        self.upsample   = nn.ConvTranspose3d(in_channels, in_channels // 2, kernel_size=2, stride=2)
        self.conv_block = ConvBlock(in_channels, out_channels)
    def forward(self, x, skip): return self.conv_block(torch.cat([self.upsample(x), skip], 1))


class nnUNet(nn.Module):
    def __init__(self, in_channels: int, out_channels: int,
                 features: Optional[List[int]] = None) -> None:
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
            self.deep_supervision_heads.append(nn.Conv3d(feat, out_channels, kernel_size=1))

    def forward(self, x):
        skips = []
        for enc in self.encoder_blocks:
            skip, x = enc(x); skips.append(skip)
        x = self.bottleneck(x)
        skips = list(reversed(skips))
        for i, dec in enumerate(self.decoder_blocks):
            x = dec(x, skips[i])
        return self.deep_supervision_heads[-1](x)


class nnUNet_VanillaUNet(nnUNet):
    def __init__(self, in_channels: int = 1, out_channels: int = 1,
                 features: Optional[List[int]] = None) -> None:
        super().__init__(in_channels, out_channels, features)


def _load_model(ckpt_path: str, features: List[int], device: torch.device) -> nnUNet_VanillaUNet:
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
# Sliding-window inference (no prior)
# ===================================================================

def infer_noprior_single_afid(
    fid: int,
    ckpt_path: str,
    img: nib.nifti1.Nifti1Image,
    patch_size: int = 64,
    batch_size: int = 7,
    overlap: float = 0.5,
    device_str: str = "cpu",
    features: Optional[List[int]] = None,
) -> Tuple[NDArray, torch.Tensor]:
    if features is None:
        features = [16, 32, 64]

    device = torch.device(device_str if (device_str == "cpu" or torch.cuda.is_available()) else "cpu")

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

    # ---- Load model ----
    t0 = time.perf_counter()
    model = _load_model(ckpt_path, features, device)
    t_load = time.perf_counter() - t0

    # ---- Gaussian map ----
    ps = patch_size
    def _g1d(n):
        s = n * 0.125; c = n // 2
        x = torch.arange(n, dtype=torch.float32)
        return torch.exp(-((x - c) ** 2) / (2 * s ** 2))
    gmap = torch.einsum("z,y,x->zyx", _g1d(ps), _g1d(ps), _g1d(ps))
    gmap = gmap / gmap.max()

    # ---- Sliding-window patch coordinates ----
    t0 = time.perf_counter()
    step = max(1, int(ps * (1.0 - overlap)))
    patch_coords = []
    for z in range(0, D, step):
        for y in range(0, H, step):
            for x in range(0, W, step):
                ze = min(z + ps, D); ye = min(y + ps, H); xe = min(x + ps, W)
                zs = max(0, ze - ps); ys = max(0, ye - ps); xs = max(0, xe - ps)
                patch_coords.append((slice(zs, ze), slice(ys, ye), slice(xs, xe)))

    num_patches = len(patch_coords)

    # ---- Extract & normalise all patches ----
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
        model_inputs.append(patch.unsqueeze(0))  # [1,D,H,W]
    t_patch = time.perf_counter() - t0

    # ---- Inference ----
    t0 = time.perf_counter()
    predictions = []
    for i in range(0, len(model_inputs), batch_size):
        chunk = torch.stack(model_inputs[i: i + batch_size]).to(device)
        with torch.inference_mode():
            preds = model(chunk)
        predictions.extend(p.cpu() for p in preds)
        del chunk, preds
    t_infer = time.perf_counter() - t0

    # ---- Gaussian-weighted reconstruction ----
    t0 = time.perf_counter()
    out = torch.zeros(vol_shape, dtype=torch.float32)
    cnt = torch.zeros(vol_shape, dtype=torch.float32)
    for pred, (zs, ys, xs) in zip(predictions, patch_coords):
        az = zs.stop - zs.start
        ay = ys.stop - ys.start
        ax = xs.stop - xs.start
        # Crop Gaussian if patch was at a boundary (< full patch_size)
        if (az, ay, ax) == (ps, ps, ps):
            w = gmap
        else:
            zo = max(0, (ps - az) // 2)
            yo = max(0, (ps - ay) // 2)
            xo = max(0, (ps - ax) // 2)
            w = gmap[zo: zo + az, yo: yo + ay, xo: xo + ax]
        pc = pred[0, :az, :ay, :ax]
        out[zs, ys, xs] += pc * w
        cnt[zs, ys, xs] += w
    out /= cnt.clamp(min=1e-8)
    t_recon = time.perf_counter() - t0

    # ---- Peak → world ----
    idx = np.unravel_index(np.argmax(out.numpy()), vol_shape)
    pred_vox = tuple(int(c) for c in idx)
    pred_world = fid_voxel2world(np.array(pred_vox, dtype=float), affine)

    t_total = t_load + t_patch + t_infer + t_recon
    print(
        f"  [noprior]  AFID {fid:02d}  patches={num_patches}  "
        f"load={t_load:.2f}s  patch={t_patch:.2f}s  "
        f"infer={t_infer:.2f}s  recon={t_recon:.2f}s  total={t_total:.2f}s  |  "
        f"pred_vox={pred_vox}  "
        f"world=[{pred_world[0]:.1f}, {pred_world[1]:.1f}, {pred_world[2]:.1f}]",
        flush=True,
    )
    return pred_world, out  # out = full-volume probability map [D,H,W]



# ===================================================================
# Snakemake entry point (All 32 AFIDs)
# ===================================================================

import os
from pathlib import Path
import csv

AFIDS_FIELDNAMES = [
    "id", "x", "y", "z", "ow", "ox", "oy", "oz", "vis", "sel", "lock",
    "label", "desc", "associatedNodeID",
]

AFID_DESCRIPTIONS = [
    "AC", "PC", "Infracollicular Sulcus", "PMJ",
    "Superior IPF", "Right Superior LMS", "Left Superior LMS",
    "Right Inferior LMS", "Left Inferior LMS", "Culmen",
    "Intermammillary Sulcus", "Right Mammilary Body", "Left Mammilary Body",
    "Pineal Gland", "Right LV at AC", "Left LV at AC",
    "Right LV at PC", "Left LV at PC", "Genu of CC", "Splenium of CC",
    "Right AL Temporal Horn", "Left AL Tempral Horn",
    "R. Sup. AM Temporal Horn", "L. Sup. AM Temporal Horn",
    "R Inf. AM Temp Horn", "L Inf. AM Temp Horn",
    "Right IG Origin", "Left IG Origin",
    "R Ventral Occipital Horn", "L Ventral Occipital Horn",
    "R Olfactory Fundus", "L Olfactory Fundus",
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
        rows.append({
            "id": lbl, "x": c[0], "y": c[1], "z": c[2],
            "ow": "0.000", "ox": "0.000", "oy": "0.000", "oz": "1.000",
            "vis": 1, "sel": 1, "lock": 1,
            "label": lbl, "desc": desc, "associatedNodeID": "",
        })
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

    pred_world, prob_map = infer_noprior_single_afid(
        fid             = fid,
        ckpt_path       = ckpt_path,
        img             = img,
        overlap = snakemake.config["inference_overlap"] if snakemake.config.get("inference_overlap") is not None else cfg_block.get("overlap", 0.5),
        patch_size      = cfg_block.get("patch_size", 64),
        batch_size      = snakemake.config.get("inference_batch_size") or 7,
        device_str      = cfg_block.get("device", "cpu"),
    )
    combined_coords[fid] = pred_world

    Path(out_coord).parent.mkdir(parents=True, exist_ok=True)
    with open(out_coord, "w") as f:
        f.write(f"{pred_world[0]} {pred_world[1]} {pred_world[2]}\n")

    Path(out_prob).parent.mkdir(parents=True, exist_ok=True)
    nib.save(
        nib.Nifti1Image(prob_map.numpy(), img.affine, img.header),
        out_prob,
    )

write_combined_fcsv(combined_coords, snakemake.output.fcsv)
