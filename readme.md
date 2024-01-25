

First make sure you have the following two packages installed:

- **devtools**
- **BiocManager**

If you don't have them yet, install them from CRAN.

Then install the **GSEAtraining** package for this course:

```r
library(devtools)
# this may take ~30 min because it builds all practice materials
install_github("jokergoo/GSEAtraining", dependencies = TRUE, 
	build_vignettes = TRUE)
```

You should also update the following two packages from GitHub because
the updates are only available on the bioc devel branch:

```r
install_github("jokergoo/simona")
install_github("jokergoo/BioCartaImage")
install_github("jokergoo/rGREAT")
```

If you have errors with installing dependency packages, try to install them manually:

```r
BiocManager::install(c("DBI", "DT", "RSQLite", "htmltools", "shiny", "AnnotationDbi", "matrixStats",
    "GO.db", "org.Hs.eg.db", "KEGGREST", "clusterProfiler", "msigdbr",
    "reactome.db", "AnnotationHub",
    "Orthology.eg.db", "microbenchmark", "ReactomePA", "DOSE", "org.Ss.eg.db",
    "CePa", "eulerr", "rGREAT", "goseq", "GSVA", "simplifyEnrichment", 
    "enrichplot", "ggplot2", "ComplexHeatmap", "circlize", "genefilter"))
install.packages("https://jokergoo.github.io/GSEAtraining_0.99.0.tar.gz", repo = NULL, type = "source")
```

The practice materials are also available at https://jokergoo.github.io/GSEAtraining/.

**Do not use it for commercial purpose.**

