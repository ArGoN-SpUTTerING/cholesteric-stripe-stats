# Data directory

This directory is intentionally empty in the git repository. The raw
measurements are distributed through Dataverse (see the DOI in the top-level
README) so that they carry a citable, versioned identifier of their own.

To run the analysis, download the two files from Dataverse and place them here:

    data/pitch_series.txt     126 x 15   manuscript Fig. 3F
    data/mixing_series.txt    135 x 12   manuscript Fig. 2D

No change to any script is needed — these are the paths the scripts expect.

## File format

Tab-separated, no header, one measurement per cell, values in micrometres.
`--` marks a cell with no measurement: only regions containing identifiable
periodic stripes were measured, so the number of measurements differs between
conditions and the columns are of unequal length.

### `pitch_series.txt` — 15 columns

Five lipid compositions, each followed by the three CB15 concentrations
(2.8 / 5 / 10 wt.-%):

    cols  1- 3   DLPC               (Fig. 3A)
    cols  4- 6   DLPC:DOPC 1:0.43   (Fig. 3B)
    cols  7- 9   DLPC:DOPC 1:1      (Fig. 3C)
    cols 10-12   DLPC:DOPC 1:2.3    (Fig. 3D)
    cols 13-15   DOPC               (Fig. 3E)

### `mixing_series.txt` — 12 columns

Three DLPC:DOPC ratios, each followed by the four total lipid concentrations
(1 / 10 / 30 / 50 uM):

    cols  1- 4   DLPC:DOPC 1:0.43
    cols  5- 8   DLPC:DOPC 1:1
    cols  9-12   DLPC:DOPC 1:2.3

The scripts do not assume this layout is correct: they verify it by checking
that within every composition the median spacing falls as the CB15
concentration rises, and that the resulting medians reproduce the boxes drawn
in the published figure.
