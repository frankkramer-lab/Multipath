if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("reactome.db")

BiocManager::install("graph")

BiocManager::install("Rgraphviz")

biopax=createBiopax(level=2)
biopax=readBiopax("pi3k.owl")
View(biopax$dt)
pi3kg=pathway2Graph(biopax,"pathway1")
