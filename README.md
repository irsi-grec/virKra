<h1 align="left">
  virKra
</h1>

**virKra** is an R package devised for importing Kraken2 output (Wood et al. 2019) to a Seurat Object from Seurat R package (Hao et al. 2023). virKra allows identifying individual cells within a Seurat Object whose DNA or RNA has been taxonomically classified by Kraken2.


## Install quick-start

virKra can be easily installed from GitHub using the following command:

```R
install.packages("remotes")
remotes::install_github("irsi-grec/virKra")
library(virKra)
```


## Usage

The following tutorial is based on the assumption that the user has already run cellranger(-arc) (Zheng et al. 2017; Satpathy et al 2019) and Kraken2 on the query FASTQ file(s), and generated a Seurat Object with single-cell information from cellranger(-arc) output.


### virKra input

For using **virKra** you need to generate for every query sample a metadata table with mandatory information. The metadata table must contain the following 7 fields (as columns):

- library_type: "Gene Expression" or "Chromatin Accessibility"
- sample: id for the current sample. Must be the same as the "orig.ident" used in Seurat.
- read1: path to R1 FASTQ file generated with Kraken2 that contains classified reads.
- read2: path to R2 FASTQ file generated with Kraken2 that contains classified reads. (optional)
- index: path to the index file from 10x.
- whitelist: path to the whitelist containing all barcodes used in the assay kit (additional information [here](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-inclusion-list-formerly-barcode-whitelist)).
- cellranger_output: path to cellranger(-arc) output for the query sample

Notice that all FASTQs from a single sample may be recapitulated in a single metadata table. Below we show an example table (metadata-virKra.csv):

| library_type  | sample | read1 | read2 | index | whitelist | cellranger_output |
| --- | --- | --- | --- | --- | --- | --- |
| Gene Expression | sample1 | /path/to/gex/kraken2/sample1_S1_L001_R2_001.fastq | | /path/to/gex/index/sample1_S1_L001_R1_001.fastq.gz | /path/to/gex/737K-arc-v1.txt.gz | /path/to/cellranger_output/sample1/outs
| Gene Expression | sample1 | /path/to/gex/kraken2/sample1_S1_L002_R2_001.fastq | | /path/to/gex/index/sample1_S1_L002_R1_001.fastq.gz | /path/to/gex/737K-arc-v1.txt.gz | /path/to/cellranger_output/sample1/outs
| Chromatin Accessibility | sample1 | /path/to/atac/kraken2/sample1_S1_L001_R1_001.fastq | /path/to/atac/kraken2/sample1_S1_L001_R3_001.fastq | /path/to/atac/index/sample1_S1_L001_R2_001.fastq.gz | /path/to/atac/737K-arc-v1.txt.gz | /path/to/cellranger_output/sample1/outs

### export classifications

First, metadata has to be loaded into the working directory:

```R
md <- loadMetadata("metadata-virKra.csv")
```

Next, FASTQ containing reads classified by Kraken2 may be loaded. For single-end reads from GEX datasets use `loadSingleEnd()` function, and for paired-end reads from ATAC datasets use `loadPairedEnd()` function. Notice that "sample" parameter may contain the same sample id as the one used in Seurat.

```R
# GEX datasets
se <- loadSingleEnd(md,
                    library_se = "Gene Expression",
                    sample = "sample1")

# ATAC datasets
pe <- loadPairedEnd(md,
                    library_pe = "Chromatin Accessibility",
                    sample = "sample1")
```

To associate reads to cells we need to recover 10x barcodes from the index file. So, barcodes are easily recovered with `recoverBarcode10x()` function:

```R
# GEX datasets
se <- recoverBarcode10x(se)

# ATAC datasets
pe <- recoverBarcode10x(pe)
```

### Output

Reads associated to cells by Kraken2 can be shown as tables using the functions `summary_counts()` and `summary_cell_counts()`. `summary_counts()` reports the number of RNA/DNA counts per cell and sample, and `summary_cell_counts()` reports the number of detected cells per sample.

```R
# GEX datasets
gex_count <- summary_counts(
  se,
  library_type = "Gene Expression")

gex_cellcount <- summary_cell_counts(
  se,
  library_type = "Gene Expression")

# ATAC datasets
atac_count <- summary_counts(
  pe,
  library_type = "Chromatin Accessibility")

atac_cellcount <- summary_cell_counts(
  pe,
  library_type = "Chromatin Accessibility")
```

`addToSeuratMetadata()` function helps to import Kraken2 classifications to a Seurat Object with single-cell information; this function accepts as input either a Seurat Object or its metadata table.

```R
# GEX datasets
pbmc <- addToSeuratMetadata(se,
                            pbmc,
                            library_type = "Gene Expression")
# ATAC datasets
pbmc <- addToSeuratMetadata(pe,
                            pbmc,
                            library_type = "Chromatin Accessibility")
```


## Disclaimer

The authors provided the information and software in good faith. Under no circumstance shall authors and IrsiCaixa have any liability for any loss or damage of any kind incurred as a result of the use of the information and software provided. The use of this tool is solely at your own risk.


## References

Wood D.E., (2019). Improved metagenomic analysis with Kraken 2. Genome Biology 20, 257. https://doi.org/10.1186/s13059-019-1891-0

Hao et al. (2023) Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology 42: 293-304. https://doi.org/10.1038/s41587-023-01767-y

Zheng G. X. Y. et al. (2017). Massively parallel digital transcriptional profiling of single cells. Nature Communications 8: 1-12. https://doi.org/10.1038/ncomms14049

Satpathy A. T. et al. (2019). Massively parallel single-cell chromatin landscapes of human immune cell development and intratumoral T cell exhaustion. Nature Biotechnology 37: 925-936. https://doi.org/10.1038/s41587-019-0206-z


## Citation

Garrido-Sanz L., et al. (2024) Optimizing detection of HIV-1 infected cells: A novel bioinformatics pipeline leveragins Kraken2 for single-cell multiomics datasets. HIV persistance during therapy (Ford Lauderdale, FL), PP1.17-54. https://doi.org/10.1016/j.jve.2024.100466

Garrido-Sanz L., et al. (2024) Optimizing detection of HIV-1 infected cells: A novel bioinformatics pipeline leveragins Kraken2 for single-cell multiomics datasets. GeSIDA (Zaragoza, Spain), P-139. https://doi.org/10.1016/S0213-005X(24)00396-3


## Acknowledgements

Thanks to all the participants, without whom this work could not be possible, and the sample processing service (IrsiCaixa). SM-L was supported by the EU Gilead Research Scholars Program in HIV award 2021, the grant 2020 BP 00046 from the Catalan Agency of Management of University and Research Grants and the grant PID2023-150316OA-I00 from the Spanish Ministry of Science, Innovation and Universities. JM-P was supported by the grant PID2019-109870RB-I00 and PID2022-139271OB-I00 from the Spanish Ministry of Science and Innovation and grant 202130-30-31-32 (NeuroCOVID) from La Marató de TV3 Foundation. Sample processing at IrsiCaixa was possible thanks to the crowdfunding initiative People in Red, and the support of CIBERINFEC, Grifols and The EU Gilead Research Scholars Program in HIV award 2021. LG-S received financial support by the grant SLT02823 000257, funded by the "Strategic plan for research and innovation in health" (PERIS), from the Catalan Department of Health. 
