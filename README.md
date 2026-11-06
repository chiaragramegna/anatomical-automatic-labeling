# AAL MNI Coordinate Label Lookup
A Python command-line script for finding the anatomical label associated with an MNI coordinate in an AAL (Automated Anatomical Labeling) NIfTI atlas.

The script reads a 3D atlas image (`.nii` or `.nii.gz`) and its matching text label file. It first checks the atlas voxel corresponding to the supplied coordinate. If that voxel has no label, it can search for the nearest labelled voxel centre within a user-defined distance.

The lookup is performed directly from the atlas files.
## Features
- Converts MNI coordinates in millimetres to voxel indices using the NIfTI affine.
- Uses zero-based voxel indexing, as required by Python and NiBabel.
- Checks the exact voxel first.
- Optionally searches for the nearest non-zero atlas voxel within a radius specified in millimetres.
- Calculates distances in world coordinates using the NIfTI affine, accounting for voxel dimensions and image orientation.
- Reads label files where the numerical code is either the first or last item on a line.
- Handles coordinates that fall outside the atlas volume and reports when no labelled voxel is found within the requested radius.
## Requirements
- Python v.3
- NumPy
- NiBabel
- SciPy
Install the dependencies with:
```bash
python -m pip install numpy nibabel scipy
```
## Input files
You need two files from the **same atlas version**:

**1. AAL NIfTI atlas**: a 3D image containing integer region codes. Background voxels should have value 0.

**2. AAL label text file**: a mapping between each integer region code and its anatomical name.

For example, the script can read label lines in either of these formats:
`1 Precentral_L
Precentral_L 1`

The label parser preserves the remaining text on the line, so labels may contain spaces. Blank lines and lines beginning with `#` or `%` are ignored.

The NIfTI image must be a discrete label atlas, not a probabilistic or statistical image. The script checks that non-zero voxel values are integer-valued.
## Configuration
Open `label_mni_coord_aal.py` and edit the default values near the top if you want to run the script without specifying all options on the command line:

```python
from pathlib import Path

AAL_PATH = Path("/path/to/aal_MNI_V4.nii")
AAL_LABELS_PATH = Path("/path/to/aal_MNI_V4.txt")

MNI_COORD = (20.0, -62.0, 54.0)
SEARCH_RADIUS_MM = 6.0
```

Replace the example paths with the locations of your own atlas and label files. You can also supply the paths and coordinate when running the script, as shown below.
## Usage
Run the script from a terminal:

```bash
python label_mni_coord_aal.py \
    --atlas "/path/to/aal_MNI_V4.nii" \
    --labels "/path/to/aal_MNI_V4.txt" \
    --coord 20 -62 54 \
    --radius 6
```

The `.nii.gz` format can also be used if supported by the selected atlas files.
## Command-line arguments

| Argument        | Required? | Description                                                                                |
| --------------- | --------- | ------------------------------------------------------------------------------------------ |
| `--atlas`       | No*       | Path to the AAL NIfTI atlas file.                                                          |
| `--labels`      | No*       | Path to the matching AAL text label file.                                                  |
| `--coord X Y Z` | No        | MNI coordinate in millimetres. The example given is `(20, -62, 54)`.                            |
| `--radius`      | No        | Maximum distance in millimetres for the nearest-labelled-voxel search. The default is `6`. |


* If omitted, the script uses the default paths configured at the top of the Python file. Those paths must be edited before running the script without `--atlas` and `--labels`.
## How the lookup works
**1.** The script loads the NIfTI atlas and its label text file.

**2.** It converts the input coordinate from world coordinates to a floating-point voxel position using the inverse NIfTI affine.

**3.** It rounds the voxel position to the nearest voxel, using MATLAB-compatible half-away-from-zero rounding.

**4.** If the voxel contains a non-zero atlas code, the matching anatomical name is printed.

**5.** If the voxel is zero or lies outside the atlas, the script finds the nearest non-zero voxel in the atlas. The distance is calculated in millimetres from the voxel centre to the original input coordinate.

**6.** If the nearest labelled voxel is no farther away than `--radius`, its code and anatomical name are printed. Otherwise, the script reports that no labelled voxel was found within the radius.

The fallback search returns the label of the **nearest labelled voxel centre**. It does not calculate the nearest region boundary, region centroid, or a majority label in a spherical neighbourhood.
## Coordinate-space note
Coordinates are interpreted as millimetres in the NIfTI image's RAS+ world-coordinate convention. The atlas must be in the intended MNI space, and the input coordinates must use the same spatial reference. The script uses the NIfTI affine; it does not register an image to MNI space or convert coordinates from another space or convention.

If the NIfTI header does not specify spatial units, the script assumes that the affine is in millimetres and issues a warning. Check the atlas metadata and coordinate convention before using the output for analysis.
## Output
The script reports:
- The input MNI coordinate.
- The corresponding voxel index (zero-based).
- Whether the exact voxel has an AAL label.
- If the fallback search is used, the nearest labelled voxel and its distance in millimetres.
- The numerical region code and its name, when the code is present in the label file.
## Limitations
- The atlas and label file must correspond to the same atlas version and coding scheme.
- The script labels one coordinate at a time.
- The nearest-voxel fallback searches all non-zero atlas voxels; it does not restrict the search to grey matter or apply tissue segmentation.
- The result depends on the atlas resolution, affine, spatial registration, and label definitions.
- A nearby label is not necessarily the anatomically correct label for a coordinate if the atlas is not correctly aligned with the coordinate space.
## Troubleshooting

**`FileNotFoundError`:**
Check the paths supplied with `--atlas` and `--labels`, or edit the default paths in the script.

**No code/name pairs could be read:**
Check that the text file contains a numerical code at the beginning or end of each label line.

**No labelled voxel found within the search radius:**
Increase `--radius` only if a larger anatomical search distance is appropriate. The output also reports the distance to the nearest labelled voxel.

**Unexpected anatomical label:**
Confirm that the coordinate and atlas use the same MNI reference and coordinate convention, and verify that the NIfTI affine and label file belong to the intended atlas version.
## Citation and attribution
If you use this script or an AAL atlas in academic work, cite the relevant atlas publication and the source from which you obtained the atlas files, following the citation guidance provided with that atlas distribution.
