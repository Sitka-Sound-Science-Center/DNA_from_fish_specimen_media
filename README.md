# An eDNA-inspired, non-destructive approach to DNA isolation from archived fish specimen fixative media

Data and R code for the statistical analysis. Everything runs from
[`Stat.Analysis.Rmd`](Stat.Analysis.Rmd); the knitted version with all output and
figures is [`Stat.Analysis.md`](Stat.Analysis.md).

## Data

`Supplementary table 1.xlsx` is the master sheet, one row per fish and 40 in all.
It has weight and length before and after fixation, tube numbers, and the
Qubit readings for both tissue and filter. The three `.txt` files are extracts
from it, and those are what the analysis reads.

`data_paired_raw.txt` is the main one: 40 fish, each extracted both ways (filter
and tissue), fixed for 24 to 672 hours at either room temperature or in the
fridge. Concentrations are ng/uL and `ID` pairs the two rows for each fish.

`data_filter.txt` is the 40 filter extractions with the fish weights attached,
lifted column for column off the master sheet. The `labels` column bins fish into
2-4, 4-6, 6-8 and 8-10 g. Its `type` column reads "tissue" on every row even
though the values are the filter measurements, and ten fish have no weight-loss
data, blank in the master sheet as well. IDs are numbered separately in each
file, so the two files do not join on `ID`.

`data_paired_norm.txt` is the tissue values normalized by fin clip weight, kept for
reference but not used.

Immersion time is written as "24h", "72h", "1 week", "2 weeks" and "4 weeks" in
the master sheet and as 24, 72, 168, 336 and 672 hours in the `.txt` files. The
lengths, tube numbers and length-loss columns exist only in the master sheet.

## Figures

`Figures/` holds the publication versions at 300 dpi PNG, vector PDF for the
single panels and TIFF for 3 and 4. `Stat.Analysis_files/figure-gfm/` holds the
copies embedded in the rendered markdown.

## Running it

```r
install.packages(c(
  "rmarkdown", "knitr", "dplyr", "ggplot2", "ggpubr", "ggpmisc", "scales",
  "lme4", "lmerTest", "pbkrtest", "emmeans", "DHARMa", "performance", "magick"
))

rmarkdown::render("Stat.Analysis.Rmd", output_format = "github_document")
```

Paths are relative to the repository root. `Stat.Analysis.md` and the
`figure-gfm` folder are both generated, so edit the `.Rmd` rather than either of
those. Run with R 4.6.1 and pandoc 3.8.3.
