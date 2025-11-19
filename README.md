# rDeamidation
This is a command-line script to calculate deamidation rates of peptides, estimating the deamidation rates of Asparagine (N) to Aspartic acid (D) and Glutamine (Q) to Glutamic acid (E) for each raw file. The program accepts an "evidence.txt" file from [MaxQuant](https://maxquant.org/) as input, and generates tab-delimited text files as output.
This repository contains an R-language port of the original Python script, [deamidation](https://github.com/dblyon/deamidation/), developed by David Lyon.

## Overview
This R script replicates the functionality of the original Python project, adapting the core logic and features to the R programming environment. Since these programs involve a stochastic process (i.e., bootstrapping), the output results are not necessarily identical.

## Usage
### Steps
1. Move the subject "evidence.txt" files (e.g., `evidence_PF01E.txt`) to a specific directory.
2. Modify file paths for the input evidence file and output directory in the `rdeamidation.R` file.
3. Execute the `rdeamidation.R` script in R environment.

### Example of input files
- `evidence_PF01E.txt`: MaxQunat "evidence.txt" file generated in an analysis of a dental enamel sample for a Pleistocene *Palaeoloxodon* specimen from Penghu Channel, Taiwan, along with an experimental blank.
- `evidence_PF02E.txt`: MaxQunat "evidence.txt" file generated in an analysis of a dental enamel sample for a Pleistocene *Sus* specimen from Penghu Channel, Taiwan, along with an experimental blank.

Please see Tsutaya et al. 2025 (Science [388:176&ndash;180](https://doi.org/10.1126/science.ads3888)) for the details of these files.

### Output files
These are essentially the same as those generated with [deamidation](https://github.com/dblyon/deamidation/). In this project, a user-defined character suffix can be added to each output file to distinguish the results from different runs.

- `Deamidation.txt`: Lists the raw files and their respective deamidation for N and Q, as mean, standard deviation, 95% confidence lower and upper limit
- `Number_of_Peptides_per_RawFile.txt`: Number of peptides available for calculation of deamidation for N and for Q
- `Bootstrapped_values.txt`: All the deamidation percentages calculated by e.g. 1000 bootstrap iterations, which are subsequently used to calculate the mean, std, and CI shown in `Deamidation.txt`
- `Protein_deamidation.txt`: Deamidation on the protein level, to be used with restraint since there are usually few data to acquire meaningful results; therefore, no bootstrapping is applied
- Deamidation on the protein level can also be calculated. However, this is NOT recommended as a default analysis, since there often is too little data for this to be meaningful/robust. Use `protein.bootstraps` argument to additionally generate the following files, which are analogous to the files above but grouped by raw files and proteins instead of only by raw files.
	- `Bootstrapped_values_proteins.txt`
	- `Deamidation_proteins.txt`
	- `Number_of_Peptides_per_Protein_per_RawFile.txt`

## Credits
- Original Python version, [deamidation](https://github.com/dblyon/deamidation/), by David Lyon
- R version, rDeamidation, by Takumi Tsutaya

## Citation
Please indicate the name of this project "rDeamidation" and cite the following original publication.
- Mackie et al. 2018. Palaeoproteomic profiling of conservation layers on a 14th Century Italian wall painting.  Angew Chem Int 57:736&ndash;7374. https://doi.org/10.1002/anie.201713020

If applicable, please also cite the released version of this script, uploaded to [Zenodo](https://zenodo.org/) with a [persistent identifier](https://doi.org/10.5281/zenodo.17595864).

## License
This project is licensed under the MIT License.  
See the [LICENSE](LICENSE) file for details.

## Acknowledgments
David Lyon agreed to publish the R version of the script. Zandra Fagern&auml;s and Bharath Nair reviewed the R version of the script.

