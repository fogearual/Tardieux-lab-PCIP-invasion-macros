# Cell-and-parasite-counting

ImageJ/Fiji macros to quantify cells and Toxoplasma gondii invasion from `.czi` images.

## Contents

Macros :
  
  - Invasion_assay_macro.ijm = Main assay. Writes a per-image QC stack, CSVs, a formatted workbook and histograms.
                           
  - Code_validation_macro.ijm = Quality control of the automated counts: randomly samples images per experiment, produces QC stacks and
                                an empty sheet to fill in manual counts, then the sheet computes the automatic-vs-manual error.

Samples : 

  - Sample_1.czi = sample czi file from our ivasion assays
  
  - result_sample_1.tiff = final QC produced by the macro

## Requirements

- macros are in µm². If your images are not calibrated, set the pixel size in Fiji first (`Image ▸ Properties`) or the area filters will be wrong.

## Input layout

The invasion pipeline detects the folder depth automatically :

- selected folder/experiment/condition/*.czi = for several experiments
- selected folder/condition/*.czi = for a single experiment
- selected folder/*czi = for a single condition

Channel colours (DIC, GFP = green, DAPI = blue, RED) are read from the CZI metadata. 
If that fails, fixed channel numbers in the macro can be modified.

## Running a macro

1. Open Fiji.
2. `Plugins ▸ Macros ▸ Run…` and pick the `.ijm` file (or drag it into the Fiji toolbar and press **Run**).
3. Select the folder when prompted. Output is written to a `Results/` folder created inside the folder you selected.

Nothing is installed and no image is modified in place — the macros only read the `.czi` files and write new files under `Results/`.

## Key parameters

The detection settings live in a clearly-marked block near the top of `Invasion_assay_macro.ijm`. 
The values below were fixed after validation and countless manual calibrations and should not be changed casually — the counting-control macro exists to check them.

## How the counting works

Each image is split into its channels, whose colours (DIC, GFP, DAPI, RED) are read from the CZI metadata. The three fluorescence channels are then segmented independently.

- Host cells (DAPI) : The blue channel is background-subtracted (rolling ball), lightly smoothed, and thresholded with the Triangle method scaled by a fixed factor; touching nuclei are separated by a watershed, and objects are kept by a size filter (≥ 105 µm²). Because nuclei often overlap in dense fields, the macro does not simply count objects: it computes the median object area for that field, and any object larger than 1.7 × that median is interpreted as a clump of several nuclei and counted as 2, 3 or 4 cells (capped at 4) depending on the area of the cluster and the median nuclei area. This makes the cell count adaptive to each field rather than to a global assumption, which matters because cell density varies from field to field.

- Green parasites (GFP). The green channel is background-subtracted, contrast-enhanced, and thresholded (Triangle × 4.25 — a deliberately stringent factor to reject the diffuse GFP background that would otherwise inflate the count), then watershed-separated. Particles are retained only if their area falls in the 7.25–62.5 µm² window ; anything smaller is noise and anything larger is debris or an unresolved cluster.

- Double-positive parasites (red and green = extracellular). The red channel is thresholded into a binary mask. For every green object, the macro measures the mean of the red mask inside that object's outline: a value near 255 means the object is fully red-positive, near 0 means not at all. A green parasite is scored as double-positive when at least 92 % of its area is red-positive.

**The quantity used for every figure is green-only parasites per 100 host cells : (green − double-positive) / nuclei × 100. Dividing by the nucleus count normalises for the cell density of each field, so a field that simply contains more cells does not appear more infected.**

**Every threshold and size window above is fixed and identical for all images and all conditions — the code applies exactly the same rules everywhere, which is the point of an automated count. It is not tuned per image. The counting validation macro exists precisely to check these settings : it re-runs this same segmentation on a random sample of images and lets you compare the automatic counts against your own manual counts on the same fields.**

## Main parameters

- `greenThresholdFactor` = 4.25 (green threshold = Triangle × factor)
- `greenMinArea` and `greenMaxArea` = 7.25 and 62.5 µm² (parasite size window)
- `redThresholdFactor` = 2.8 (red threshold = Triangle × factor)
- `minOverlapPercent` = 92 (% of a green object that must be red-positive to count as double-positive)
- `nucleiDoubletFactor` = 1.7 (object larger than this × median is considered as cluster)

**Every other modifiable parameter is in section 2 of the macro : "detection"**

## Output

- `Results/Stacks_<name>/` = one 5-slice QC stack per analyzed image (blue, green,  red, red+green merge, full merge), with outlines and counts burnt in.
- `Results/Data/` = Raw data for Image J : nucleus ROIs and one CSV per analyzed condition.
- `Results/Final results <name>.xls` = per-experiment summary (top) and per-image detail (bottom).
- `Results/Histogram <name>.png` for every experiment and `Histogram pooled <name>.png` for the total.

## A note on statistics

The 50 fields recorded per condition are pseudo-replicates, they come from one coverslip. 
The independent unit is the experiment. Error bars in the pooled (total) figures are therefore computed with the experiment as the unit (n = number of experiments), not the field. 

## Code valdation macro

The quality control macro randomly picks czi files, produces QC stacks with the exact same process as the invasion assay macro and creates an empty sheet to fill in manual counts, then the sheet computes the automatic-vs-manual error. The stack contains the raw blue channel for manual cell counting, then the raw green channel, raw red channel, and raw merged green + red channels for manual parasite counting, then the auto-counting annotated channels for validation.

## Citation

If you use this code, please cite it — see `CITATION.cff` file.

## License

MIT — see `LICENSE`.
