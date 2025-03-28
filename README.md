# Sleeping Beauty Paper

This repository contains all analysis code and data for the paper "A forward genetic screen identifies Sirtuin1 as a driver of neuroendocrine prostate cancer."

## TAPDANCE Analysis (Perl)

The TAPDANCE (Transposon Analysis Pipeline for Driver And Cancer Element identification) implementation for this project is in the `TAPDANCE_NEPC` subdirectory.

Full instructions for the TAPDANCE pipeline, including setup, configuration, and execution, can be found in the [TAPDANCE README](TAPDANCE_NEPC/README.md).

The large input and intermediate files not included in this repository are available at [Zenodo](https://zenodo.org/records/15098998). Final common insertion site level data outputs are contained in the `TAPDANCE_NEPC` subdirectory.

## Downstream Analysis (R)

The main.sh file contains a bash pipeline that will recreate the entire analysis.

To run this pipeline on a Mac, you can do the following:

```bash
bash PATH_TO_DIR/sleeping-beauty-paper/main.sh
```

This will create an `experiments/` directory within `sleeping-beauty-paper/` containing the results.

### Requirements

Note: pandoc is required to generate the network figures. When running the code within RStudio, it is loaded through the application. However, when running it with an RScript (which is what our bash pipeline does), you may encounter the following error:

```
Error in htmlwidgets::saveWidget(graph, file, selfcontained = selfcontained, :
    Saving a widget with selfcontained = TRUE requires pandoc.
```

To address this, include pandoc in your environment's PATH variable. In Mac, this can be done by adding the following lines of code into your shell source (e.g. .zshrc) as shown below:

```bash
PATH="${PATH}:/Applications/RStudio.app/Contents/MacOS/quarto/bin/tools/pandoc"
export PATH
```

After reloading your shell, you should now be able to run the script without encountering any errors.
