# sleeping-beauty-paper

The main.sh file contains a bash pipeline that will recreate the entire analysis.

To run this pipeline on a Mac, one can do the following:

```
bash PATH_TO_DIR/sleeping-beauty-paper/main.sh
```

This will create a `experiments/` directory within `sleeping-beauty-paper/` containing the results.

Code to generate figures for the Sleeping Beauty Paper

Note: pandoc is required to generate the network figures. When running the code within RStudio, it is loaded through the application. However, when running it with an RScript (which is what our bash pipeline does), you may encounter the following error:

Error in htmlwidgets::saveWidget(graph, file, selfcontained = selfcontained,  :
    Saving a widget with selfcontained = TRUE requires pandoc.

To address this, include pandoc in your environment's PATH variable. In Mac, this can be done by adding the following lines of code into your shell source (e.g. .zshrc) as shown below:

PATH="${PATH}:/Applications/RStudio.app/Contents/MacOS/quarto/bin/tools/pandoc"
export PATH

After reloading your shell, you should now be able to run the script without encountering any errors.

