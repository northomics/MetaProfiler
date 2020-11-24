# Contents of this file

 * Introduction
 * Citation
 * Installation
 * Examples
 * License

# Introduction
![alt text](https://github.com/psmyth94/MetaProfiler/blob/master/man/logo/logo.png)

This repository contains the source code of the `MetaProfiler` software package.
It provides calculations for local false discovery rates of protein-based stable isotopic probing (SIP) results and performs taxonomic, functional, phylogenetic, and time series analysis of microbiome dynamics.

`MetaProfiler` has only been tested on MetaProSIP results from OpenMS, but it is designed to work with multiple tools that extract heavy peptide features from light peptide identifications. More tools such as ProteinTurnover will be tested in the future.

# Citation

If you use `MetaProfiler` in your projects, please cite the preprint

Patrick Smyth, Xu Zhang, Zhibin Ning, Janice Mayne, Jasmine I Moore, Krystal Walker, Mathieu Lavallée-Adam and Daniel Figeys 2020, *Studying the dynamics of the gut microbiota using metabolically stable isotopic labeling and metaproteomics* [doi:10.1101/982884](https://doi.org/10.1101/2020.03.09.982884)

# Installation

Make sure to have `R >= 3.5.0` installed. Paste the following lines into your `R` session.

```{R}
# instal devtools, if you do not have it.
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# install MetaProfiler via devtools
library(devtools)
install_github("psmyth94/MetaProfiler")
```

# Examples

This an example R script that creates the MetaProfiler class object.

```{R}
  
# the units for the time measurements.
time_unit = "hour"
# the start time of when the microbiome was fed metabolic stable isotopes.
time_zero = 0
# name of the incorporation measurements. MetaProSIP calls it RIA.
incorporation_name = "RIA"
# name of the intensity measurements. MetaProSIP calls it INT.
intensity_name = "INT"
# name of the score values. MetaProSIP uses correlation scores.
score_name = "Cor."
# name of the labeling ratio measurements. We will not specify it yet as we will not use the LR values from MetaProSIP.
labeling_ratio_name = NULL
# automatically create the experimental design table using the names in the result directory.
design <- create_experimental_design(
  results_directory = "./Protein_SIP_results", # the directory with the files containing the information about the heavy peptide features.
  Sample = "Ref\\d+", # the sample names
  hour = "(?<=_D)\\d+|(?<=_D\\-)\\d+" # the time when the sample was collected. Be sure to use the same name as specified in variable time_unit.
)

# Create MetaProfiler class object.
Object <- MetaProfiler(
  design,
  # The file names containing the information about the heavy peptide features.
  # Does not need to be specified if the `design` contain a column with the filenames.
  data = NULL,
  time_unit = time_unit,
  time_zero = 0,
  # the names of the variables.
  incorporation_name = "RIA",
  intensity_name = "INT",
  labeling_ratio_name = NULL,
  labeling_ratio_columns = NULL,
  score_name = "Cor.",
  # the columns containing the information.
  # Does not need to be specified as long as the first word is the corresponding variable's name, followed by a unique identifier
  # e.g. [RIA 1, RIA 2, RIA 3, etc] or [RIA light, RIA heavy].
  accession_column = "Protein Accessions"
)

# all files were generated using MetaLab (http://dashboard.imetalab.ca/#/).
# The pep2pro, pep2taxon, and pro2func variable can also be a matrix, data.frame, or data.table.
# The function will try to guess the accession, peptide, taxon, and function columns of these tables.
# The peptide and taxon columns are relatively easy to guess, but the accession and function columns can be a little tricky.
# For the accession column, by default, `make_annotation_table` will look for columns with uniprot accession patterns.
# For the function column, it will look for column names containing typical functional annotation databases such as KEGG, BRITE, GO, COG, and NOG.
# If it fails to guess, you will need to specify the column names.
Object <- make_annotation_table(
  Object,
  pep2pro = pep2pro, # the peptide to protein file.
  pro2func = pro2func, # the accession to function file.
  pep2taxon = pep2taxon, # the peptide to taxon file.
  compute_razor_protein = T # compute razor protein. See https://med.uottawa.ca/core-facilities/facilities/proteomic/resources/interpret-result for details about what are razor proteins.
)

# Extract the heavy peptide features
Object <- get_data(
  Object,
  light_peptide = F # Since the samples do not have a light protein spike-in. We will not assume that the light peptide will always be present.
)

# Compute an LFDR for each heavy peptide features
Object <- lfdr(Object)

# Filter at an LFDR of 1/(1 + (1/10)) = 0.0909. This means that the odds for false discovery is 1 out of 10.
Object <- filter_data(Object, LFDR_threshold = "Strong")

# Check the distribution of LR.
plot(Object, "LR")

# Check the distribution of RIA.
plot(Object, "RIA")										  
```

# Licensing and contributions
`MetaProfiler` is licensed under the GPL (>= 2) license. Contributions are welcome.