# Multipath
Creating integrated reproducible pathway multilayer models

![alt text](https://github.com/frankkramer-lab/Multipath/blob/master/img/multipath.png "multipath")
## Introduction
Biological pathway data integration has become a topic of interest in the past years. This interest originates essentially from the continuously increasing size of existing prior knowledge as well as from the many challenges scientists face when studying biological pathways. Multipath is a framework that aims at helping re-trace the use of specific pathway knowledge in specific publications, and easing the data integration of multiple pathway types and further influencing knowledge sources.
Using Multipath, BioPax-encoded pathways can be parsed and embedded into multilayered graphs. Modifications can be applied to these graphs to generate different views.
The package is implemented as a part of [the Multipath Project](https://www.sys-med.de/en/junior-research-groups/multipath/)  directed by [Dr. Frank Kramer](https://www.uni-augsburg.de/de/fakultaet/fai/informatik/prof/misit/mitarbeiter/) .
## Installation
### Preinstallation
Multipath depends on multiple packages. The packages are the following:
[UniProt.ws](https://www.bioconductor.org/packages/release/bioc/html/UniProt.ws.html), [dbparser](https://cran.r-project.org/web/packages/dbparser/index.html), [rBiopaxParser](https://github.com/frankkramer-lab/rBiopaxParser), [mully](https://github.com/frankkramer-lab/mully), [TCGAretriever](https://cran.r-project.org/web/packages/TCGAretriever/index.html), [stringr](https://cran.r-project.org/web/packages/stringr/), [svMisc](https://cran.r-project.org/web/packages/svMisc/), [uuid](https://cran.r-project.org/web/packages/uuid/), [dplyr](https://cran.r-project.org/web/packages/dplyr/), [crayon](https://cran.r-project.org/web/packages/crayon/)
Please make sure to install the packages [UniProt.ws](https://www.bioconductor.org/packages/release/bioc/html/UniProt.ws.html), [rBiopaxParser](https://github.com/frankkramer-lab/rBiopaxParser) and [mully](https://github.com/frankkramer-lab/mully) before using the package.

To install the [UniProt.ws](https://www.bioconductor.org/packages/release/bioc/html/UniProt.ws.html) package:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("UniProt.ws")
```
To install the [mully](https://github.com/frankkramer-lab/mully) package:
```R
require(devtools)
install_github("frankkramer-lab/mully")
library(mully)
```
To install the [rBiopaxParser](https://github.com/frankkramer-lab/rBiopaxParser) package:
```R
require(devtools)
install_github("frankkramer-lab/rBiopaxParser")
library(rBiopaxParser)
```

### Installation via Github

```R
require(devtools)
install_github("frankkramer-lab/Multipath")
library(Multipath)
```
## Test the package
In this section, we provide a demo to test the package by calling some of the function. To run the script, you need to download the [Signaling by Wnt from the Reactiome database](https://reactome.org/content/detail/R-HSA-195721) in the BioPax format.
### Create a mully graph from a BioPax-encoded pathway
```R
wntBiopax=readBiopax("wntpathway_reactome.owl")
pathwayID=listPathways(wntBiopax)$id[1]
wntmully=pathway2Mully(wntBiopax,pathwayID)
plot3d(wntmully,layers=T,vertex.label=NA,edge.width=5,edge.arrow.size=5)
```
![alt text](https://github.com/frankkramer-lab/Multipath/blob/master/img/wntpathway.png "Wnt pathway mully model")

### Generate Views
```R
view1=pathwayView(wntmully,"View1")
view1=addStep(view1,action = "remove",element="layer",name="Rna",trans = T)
suppressWarnings(plot3d(view1$modified))
```
![alt text](https://github.com/frankkramer-lab/Multipath/blob/master/img/view1rna.png "view1")
```R
view2=pathwayView(wntmully,"View2")
view2=addStep(view2,action = "remove",element="layer",name="PhysicalEntity",trans=T)
suppressWarnings(plot3d(view2$modified))
```
![alt text](https://github.com/frankkramer-lab/Multipath/blob/master/img/view2physentity.png "view2")
```R
view3=pathwayView(wntmully,"View3")
view3=addStep(view3,action = "remove",element="layer",name="Complex",trans=T)
suppressWarnings(plot3d(view3$modified))
```
![alt text](https://github.com/frankkramer-lab/Multipath/blob/master/img/view3complex.png "view3")

## Available Functions
Multipath functions are divided into different files depending on their functionnality range:
[Reactome](https://github.com/frankkramer-lab/Multipath/blob/master/R/Reactome.R) ,
[Views' Functions](https://github.com/frankkramer-lab/Multipath/blob/master/R/Views.R) ,
[DrugBank Functions](https://github.com/frankkramer-lab/Multipath/blob/master/R/DrugBank.R) ,
[UniProt Functions](https://github.com/frankkramer-lab/Multipath/blob/master/R/UniProtKB.R) ,
[Intraedges between nodes Functions](https://github.com/frankkramer-lab/Multipath/blob/master/R/DrugBank.R) ,
[Demo Functions](https://github.com/frankkramer-lab/Multipath/blob/master/R/wnt_pathway.R) ,
[Integrated Model Function](https://github.com/frankkramer-lab/Multipath/blob/master/R/Multipath.R)


| Function |Description|
| --------------- |-----------|
|`downloadPathway(pathwayID,biopaxLevel,destDirectory,overwrite)`|Download a BioPax encoded Reactome pathway function|
|`pathway2Mully(biopax,pathwayID)`|Build a mully graph from a BioPax encoded pathway function|
|`pathwayView(g,name)`|Constructor Function, Create an empty view|
|`print(v)`|Print function|
|`addStep(v,action,element,name,layerName,V1,V2,attributesnames,attributes,multi,trans)`|Document the modification of a graph function|
|`undo(v,stps)`|Undo stps number of steps from the view v|
|`loadDBXML(file)`|Parse DrugBank XML file function|
|`getDBDrug(data,drug)`|Get a DrugBank entry from DrugBank function|
|`getDBDrugInteractions(data,drug)`|Get Drugs interactions from DrugBank function|
|`getDBDrugEnzymes(data,drugList)`|Get Drug-Enzymes relations from DrugBank function|
|`getDBDrugTransporters(data,drugList)`|Get Drug-Transporters relations from DrugBank function|
|`getDBDrugCarriers(data,drugList)`|Get Drug-Carriers relations from DrugBank function|
|`getDBDrugTargets(data,drugList)`|Get Drug-Targets relations from DrugBank function|
|`addDBLayer(g,data,drugList)`|Add a DrugBank Layer to a graph function|
|`getDBtoUPKB(data,drugLiy,proteinList)`|Get Drug-Proteins relations from DrugBank function|
|`getAllUPKB(up)`|Get all UniProt protein entries function|
|`getUPKBInfo(up,proteins,col)`|Get UniProt protein entry infos function|
|`getUPKBInteractions(up,proteins)`|Get Proteins interactions from UniProt function|
|`getUPKBtoDB(up,proteinList,drugList)`|Get Proteins-Drug relations from UniProt function|
|`addUPKBLayer(g,up,proteinList)`|Add a UniProt Layer to a graph function|
|`getUPKBDBRelations(up,data,proteinList,drugList)`|Get Proteins-Drug relations from UniProt and DrugBank function|
|`multipath(name,up,proteinList,data,drugList)`|Build integrated model function|
|`wntpathway(file)`|Track and undo demo function|
