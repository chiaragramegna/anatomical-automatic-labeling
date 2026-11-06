#!/usr/bin/env python3
"""Look up an AAL anatomical label for an MNI coordinate.

Dependencies:
python -m pip install numpy nibabel scipy

NIfTI voxel indices in Python are zero-based. MNI coordinates are supplied
in millimetres in the NIfTI image's RAS+ world-coordinate convention.
"""

from __future__ import annotations

import argparse
import math
import warnings
from pathlib import Path
import typing

import numpy as np # type: ignore

# ---------------------------------------------------------------------------------
# Defaults: edit these paths according to your local installation of the AAL atlas.
# ---------------------------------------------------------------------------------
AAL_PATH = Path("/.../aal_MNI_V4.nii")
AAL_LABELS_PATH = Path("/.../aal_MNI_V4.txt")

#---------------------------------------------------------------------------------
# Example of MNI coordinate and search radius (mm) for nearest-region fallback.
# Edit as needed, or supply your own values on the command line.
#---------------------------------------------------------------------------------
MNI_COORD = (20.0, -62.0, 54.0)
SEARCH_RADIUS_MM = 6.0


def _integer_code(token: str):
    """Return an integer for an integer-valued numeric token; otherwise None."""
    try:
        value = float(token)
    except ValueError:
        return None

    if math.isfinite(value) and value.is_integer():
        return int(value)
    return None


def load_labels(labels_path: Path) -> typing.Dict[int, str]:
    """Read code/name pairs, allowing the numeric code at either end of a line.

    Supported examples:
        1 Precentral_L
        Precentral_L 1

    The complete label text is retained, including embedded spaces.
    Blank lines and lines beginning with # or % are ignored.
    """
    if not labels_path.is_file():
        raise FileNotFoundError(f"Label file not found: {labels_path}")

    labels: typing.Dict[int, str] = {}
    with labels_path.open("r", encoding="utf-8-sig", errors="replace") as label_file:
        for line_number, raw_line in enumerate(label_file, start=1):
            line = raw_line.split("#", 1)[0].strip()
            if not line or line.startswith("%"):
                continue

            fields = line.split()
            if len(fields) < 2:
                continue

            first_code = _integer_code(fields[0])
            last_code = _integer_code(fields[-1])

            if first_code is not None:
                code = first_code
                label = " ".join(fields[1:])
            elif last_code is not None:
                code = last_code
                label = " ".join(fields[:-1])
            else:
                # Ignore headers and other lines that do not contain a code
                # in the first or last column.
                continue

            if label:
                labels[code] = label

    if not labels:
        raise ValueError(
            f"No code/name pairs could be read from {labels_path}. "
            "Check the text-file format."
        )

    return labels


def round_half_away_from_zero(values: np.ndarray) -> np.ndarray:
    """Use round() behavior, including exact half-way values."""
    values = np.asarray(values, dtype=float)
    rounded = np.where(values >= 0, np.floor(values + 0.5), np.ceil(values - 0.5))
    return rounded.astype(np.int64)


def _affine_in_mm(image) -> np.ndarray:
    """Return the image affine with its world coordinates expressed in mm."""
    affine = np.asarray(image.affine, dtype=float).copy()
    if affine.shape != (4, 4) or not np.all(np.isfinite(affine)):
        raise ValueError("The NIfTI image does not have a valid 4x4 affine matrix.")

    spatial_unit = image.header.get_xyzt_units()[0]
    unit_to_mm = {
        "meter": 1000.0,
        "mm": 1.0,
        "micron": 0.001,
        None: 1.0,
        "unknown": 1.0,
    }
    if spatial_unit not in unit_to_mm:
        raise ValueError(f"Unsupported NIfTI spatial unit: {spatial_unit!r}")
    if spatial_unit in (None, "unknown"):
        warnings.warn(
            "NIfTI spatial units are unknown; assuming the affine is in mm.",
            RuntimeWarning,
            stacklevel=2,
        )

    # Scale all three world-coordinate rows, including their translations.
    affine[:3, :] *= unit_to_mm[spatial_unit]

    if abs(np.linalg.det(affine[:3, :3])) < 1e-12:
        raise ValueError("The NIfTI affine is singular and cannot be inverted.")

    return affine


def _as_label_code(value: float) -> int:
    """Convert an atlas voxel value to an integer label code."""
    rounded = int(np.rint(float(value)))
    if not np.isclose(float(value), rounded, atol=1e-4):
        raise ValueError(
            f"Atlas contains a non-integer label value ({value}). "
            "Make sure this is a discrete label atlas, not a probability map."
        )
    return rounded


def lookup_aal_label(
    atlas_path: Path,
    labels_path: Path,
    mni_coord: typing.Iterable[float],
    search_radius_mm: float,
) -> None:
    """Print the label at an MNI coordinate, with a nearest-region fallback."""
    try:
        import nibabel as nib
    except ImportError as exc:
        raise ImportError(
            "nibabel is required. Install dependencies with: "
            "python -m pip install numpy nibabel scipy"
        ) from exc

    try:
        from scipy.spatial import cKDTree
    except ImportError as exc:
        raise ImportError(
            "SciPy is required for the millimetre-accurate nearest-region search. "
            "Install it with: python -m pip install scipy"
        ) from exc

    atlas_path = Path(atlas_path)
    labels_path = Path(labels_path)
    coord = np.asarray(tuple(mni_coord), dtype=float)

    if coord.shape != (3,) or not np.all(np.isfinite(coord)):
        raise ValueError("MNI coordinate must contain three finite values: x, y, z.")
    if not math.isfinite(search_radius_mm) or search_radius_mm < 0:
        raise ValueError("Search radius must be a finite, non-negative number of mm.")
    if not atlas_path.is_file():
        raise FileNotFoundError(f"Atlas file not found: {atlas_path}")

    labels = load_labels(labels_path)

    image = nib.load(str(atlas_path))
    data = np.asanyarray(image.dataobj)
    if data.ndim != 3:
        raise ValueError(f"Expected a 3D atlas, but the image has shape {data.shape}.")

    affine_mm = _affine_in_mm(image)

    # Nibabel uses zero-based voxel indices. The affine maps those indices to
    # world coordinates in mm; solve the inverse mapping for the input point.
    voxel_float = np.linalg.solve(affine_mm, np.append(coord, 1.0))[:3]
    voxel_index = round_half_away_from_zero(voxel_float)

    print(
        f"MNI coordinate [ {coord[0]:g}, {coord[1]:g}, {coord[2]:g} ] mm "
        f"-> voxel [ {voxel_index[0]}, {voxel_index[1]}, {voxel_index[2]} ] (0-based)"
    )

    shape = np.asarray(data.shape, dtype=np.int64)
    inside_volume = bool(np.all(voxel_index >= 0) and np.all(voxel_index < shape))

    if inside_volume:
        exact_value = float(data[tuple(voxel_index)])
    else:
        exact_value = 0.0
        print("The rounded voxel is outside the atlas volume.")

    if exact_value != 0.0:
        code = _as_label_code(exact_value)
        print("The requested voxel belongs to an AAL region.")
    else:
        print(
            "The requested voxel has no AAL label (or is outside the atlas). "
            f"Searching for the nearest labelled voxel within {search_radius_mm:g} mm..."
        )

        # Build a spatial index from the centres of every non-zero atlas voxel.
        # Distances are evaluated in MNI/world mm, so voxel size and orientation
        # are handled by the NIfTI affine rather than approximated in voxel space.
        nonzero_mask = np.isfinite(data) & (data != 0)
        voxel_coords = np.argwhere(nonzero_mask)
        if voxel_coords.size == 0:
            print("The atlas contains no non-zero labelled voxels.")
            return

        linear = affine_mm[:3, :3]
        translation = affine_mm[:3, 3]
        world_coords_mm = voxel_coords @ linear.T + translation

        tree = cKDTree(world_coords_mm)
        distance_mm, nearest_index = tree.query(coord, k=1)

        if float(distance_mm) > search_radius_mm:
            print(
                f"No labelled voxel found within {search_radius_mm:g} mm. "
                f"The nearest labelled voxel is {float(distance_mm):.2f} mm away."
            )
            return

        nearest_voxel = voxel_coords[int(nearest_index)]
        code = _as_label_code(data[tuple(nearest_voxel)])
        nearest_mni = world_coords_mm[int(nearest_index)]

        print(
            "Nearest labelled voxel: "
            f"[ {nearest_voxel[0]}, {nearest_voxel[1]}, {nearest_voxel[2]} ] (0-based)"
        )
        print(
            "Nearest voxel centre in MNI mm: "
            f"[ {nearest_mni[0]:.2f}, {nearest_mni[1]:.2f}, {nearest_mni[2]:.2f} ]"
        )
        print(f"Distance from requested coordinate: {float(distance_mm):.2f} mm")

    region_name = labels.get(code)
    if region_name is not None:
        print(f"Region code: {code} -> {region_name}")
    else:
        print(f"Region code: {code} (name not found in the label file)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Find the AAL anatomical label for an MNI coordinate."
    )
    parser.add_argument("--atlas", type=Path, default=AAL_PATH, help="AAL NIfTI file")
    parser.add_argument("--labels", type=Path, default=AAL_LABELS_PATH, help="AAL text label file")
    parser.add_argument(
        "--coord",
        type=float,
        nargs=3,
        metavar=("X", "Y", "Z"),
        default=MNI_COORD,
        help="MNI coordinate in mm (default: %(default)s)",
    )
    parser.add_argument(
        "--radius",
        type=float,
        default=SEARCH_RADIUS_MM,
        help="Maximum fallback search distance in mm (default: %(default)s)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        lookup_aal_label(args.atlas, args.labels, args.coord, args.radius)
    except (FileNotFoundError, ImportError, ValueError, OSError) as exc:
        raise SystemExit(f"Error: {exc}") from exc
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
