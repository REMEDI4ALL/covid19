# Covid: Drug Combinations

This folder contains Cell Painting results for two-drug combinations against COVID-19, including some analysis code.

Data consist of the [CellProfiler](https://github.com/CellProfiler/CellProfiler.git) features that were first extracted for each cell and subsequently were averaged (by median) within each well. Among the resulting averaged features, we have four main types of measurements (cells):
1. `combo`: cells treated with a two-drug combination;
2. `single`: cells treated with a single drug: antiviral control (GC376, Nirmatrelvir, Remdesivir) or positive control (Fenbendazole, Fluphenazine Dihydrochloride, Etoposide, Staurosporine);
2. `infected`: cells exposed to Covid but not treated with any drug;
3. `uninfected`: healthy cells.

## Experiment
- number of unique drug compounds: 184;
- number of unique drug combinations: 924;
- cell line: A549;
- solvent: Dimethyl Sulfoxide (DMSO);
- plate type: PE-6057302;
- time seeded: 2024-09-23 09:00:00;
- time painted: 2024-09-26 09:00:00;
- cells per well: ~4000.

Each drug compound was tested at two concentrations. Therefore, we have 2\*2=4 measurements corresponding to the same pair of drug compounds, each having multiple replicates for additional validation.


## Files
- `data/covid_combo.parquet.gzip`: data file (table), where each row represents morphological features obtained from Cell Painting -> CellProfiler features, per cell -> Median-averaged features, per well. First five columns contain metadata. Can be read with, e.g., `pandas.read_parquet()`;
- `data/id_to_name.txt`: two-column correspondence {batch_id: compound name}. We do not use compound names in the data file;
- `Covid Combo.ipynb`: python notebook with analysis based on Mahalanobis distance, i.e., ranking of combinations based on the similarity between morphological profile of a cell treatead by some drug combination and uninfected cells;
- `mahalanobis.py`: functions and utilities for `Covid Combo.ipynb`;
- `requirements.txt`: python requirements for analysis code.


## Metadata columns
In addition to columns corresponding to morphological features from CellProfiler, we have metadata columns:
1. `type`: {'single', 'combo', 'infected', 'uninfected'}, i.e., single drug (control), drug combination, DMSO, and healthy;
2. `batch_id_1`: unique identifier of the first drug in the given combination (find the corresponding compound name in the file `data/id_to_name.txt`);
3. `batch_id_2`: id of the second drug. If not a combination, then `batch_id_2`=`batch_id_1`;
4. `conc_1`: concentration of the first drug in the given combination, in micromoles;
5. `conc_2`: concentration of the second drug. If not a combination, then `conc_2`=`conc_1`.

**Note**: although we use the indexes `_1` and `_2` to distinguish two drugs in a combination, it does not imply any order –– drug combinations are 'symmetric' in this sense.