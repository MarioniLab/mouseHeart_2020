### Single-cell RNA-seq atlas of the earliest stages of mouse heart development

Code used to analyse the single-cell RNA-seq data from *Characterization of a common progenitor pool of the epicardium and myocardium*, Tyser, Ibarra-Soria et al., *Science*, 2021 [DOI: 10.1126/science.abb2986](https://doi.org/10.1126/science.abb2986).

The repository contains:

- `scripts`: contains all `R markdown` scripts used to analyse the data and produce figures. Output in `md` and `html` is provided to explore the results.

- `data`:  contains a file with gene annotation for the mouse reference genome, based on Ensembl version 87. To run the scripts, it is necessary to download and add to this folder the count matrix and metadata for the unbiased and reference datasets. Both are provided as supplementary information with the paper.
    + `DataS1.csv`: contains the raw counts for all `unbiased` samples. Download from [here](https://science.sciencemag.org/highwire/filestream/755147/field_highwire_adjunct_files/0/abb2986_DataS1.csv). Read this file into R and save as an `Rds` file, named `heartData_unbiased.RAW.Rds`. Alternatively, download the R object directly from [here](https://content.cruk.cam.ac.uk/jmlab/mouseEmbryonicHeartAtlas/heartData_unbiased.RAW.Rds).
    + `DataS3.csv`: contains the metadata for all `unbiased` samples. Download from [here](https://science.sciencemag.org/highwire/filestream/755147/field_highwire_adjunct_files/2/abb2986_DataS3.csv).
    + `DataS2.csv`: contains the raw counts for all `reference` samples. Download from [here](https://science.sciencemag.org/highwire/filestream/755147/field_highwire_adjunct_files/1/abb2986_DataS2.csv). Read this file into R and save as an `Rds` file, named `mesodermData_reference.RAW.Rds`. Alternatively, download the R object directly from [here](https://content.cruk.cam.ac.uk/jmlab/mouseEmbryonicHeartAtlas/mesodermData_reference.RAW.Rds).
    + `DataS4.csv`: contains the metadata for all `reference` samples. Download from [here](https://science.sciencemag.org/highwire/filestream/755147/field_highwire_adjunct_files/3/abb2986_DataS4.csv).

- `shinyApp`: scripts used for the [shiny app](https://marionilab.cruk.cam.ac.uk/heartAtlas/).
