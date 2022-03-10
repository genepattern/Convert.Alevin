<br>
<img src="GP_logo.png" alt="GenePattern logo" width="300"/>

# Convert.Alevin Documentation

**Description:** Read alevin output files and convert them to anndata (.h5ad) or loom (.loom) format
for downstream analysis. Designed to convert alevin files from runs with a velocity processed
transcriptome into a compatible format for scVelo.

**Author:** Anthony S. Castanza \
**Contact:** [genepattern.org/help](https://genepattern.org/help)

**Summary:** In order to use alevin quantifications from the alevin velocity pipeline with downstream
tools like scVelo, the quantifications must be converted into anndata or loom format containing
separate layers for the spliced and unspliced matrices. This module accepts alevin quantification
directories (in tar.gz format) and produces scVelo compatible anndata or loom files.

**Parameters:**

| Name           | Description                                                                                                                                                                                                                                                                                                 |
|----------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Alevin Quants* | Alevin output directories. Must be in .tar.gz format, and have unique names. The directory name will be used to annotate the cell barcode origins in the resulting output files.                                                                                                                            |
| Features       | A tab delimited file containing the list of spliced gene ids in column 1, the unspliced gene ids in column 2, and gene names (symbols) in column 3. (From PreprocessVelocityTranscriptome, or manually created). Optional. If not provided the module will assume unspliced genes end with the suffix “-I”. |
| Merge*         | Merge (outer join) multiple alevin quantifications into a single output file (true) or produce separate output files for each input file (false). If “true” the input directory name will be appended to the cell barcodes and stored as a “batch” layer for downstream analysis. <br> Default = true       |
| Out Type*      | Output anndata h5ad files (h5ad) or loom files (loom) <br> Default = h5ad                                                                                                                                                                                                                                   |

**Output File(s):** one or more h5ad or loom files

**Module Language:** Python \
**Source Repository:** https://github.com/genepattern/Convert.Alevin/releases/tag/v1 \
**Docker image:** [genepattern/convert_alevin:v1](https://hub.docker.com/layers/196702866/genepattern/convert_alevin/v1/images/sha256-9eee01313b073a752fcbfa4f80b67f6a7e844949610e4dbac2e809b188b19e5c?context=repo)

| Version | Comment          |
|---------|------------------|
| 1       | Initial release. |