### vignettes

This folder contains a collection of vignettes for the online R package documentation. They are not part of the R package itself. To add or update an entry:

1. Modify/add the `.Rmd` R Markdown file.
2. Update the vignette entry layout described under "Tutorials" in `_pkgdown.yml`. This list is alphabetical.
4. Make sure that it runs: render with R Markdown or use `pkgdown::build_articles`. Plots should be limited to one per `{r}` code segment in order to render properly.
5. The vignettes are built using the online continuous integration, they are not tested on a Pull Request. If any additional packages are needed for the vignettes (ideally **not**), they should be installed here, see the continuous integration setup for examples.
