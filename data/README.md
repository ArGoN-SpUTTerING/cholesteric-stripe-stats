# Data directory

This directory is intentionally empty in the git repository. The raw
measurements are distributed through Dataverse (see the DOI in the top-level
README) so that they carry a citable, versioned identifier of their own.

To run the analysis, download this file from Dataverse and place it here:

    data/pitch_series.txt     126 x 15   manuscript Fig. 3F

No change to any script is needed — this is the path the scripts expect.

## File format

Tab-separated, no header, one measurement per cell, values in micrometres.
`--` marks a cell with no measurement: only regions containing identifiable
periodic stripes were measured, so the number of measurements differs between
conditions and the columns are of unequal length.

Five lipid compositions, each followed by the three CB15 concentrations
(2.8 / 5 / 10 wt.-%):

    cols  1- 3   DLPC               (Fig. 3A)
    cols  4- 6   DLPC:DOPC 1:0.43   (Fig. 3B)
    cols  7- 9   DLPC:DOPC 1:1      (Fig. 3C)
    cols 10-12   DLPC:DOPC 1:2.3    (Fig. 3D)
    cols 13-15   DOPC               (Fig. 3E)

The scripts do not assume this layout is correct: they check that within every composition the
median spacing falls as the CB15 concentration rises.
