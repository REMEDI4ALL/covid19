# Covid: Drug Combinations

This folder contains Cell Painting results for drug combinations against COVID-19, including some analysis code.

## Experiment
Configurations:
- \# of drug combinations: 924;
- \# of drug compounds: 184;
- cell line: A549;
- solvent: Dimethyl Sulfoxide (DMSO);
- plate type: PE-6057302;
- time seeded: 2024-09-23 09:00:00;
- time painted: 2024-09-26 09:00:00;
- cells per well: 4000.

Each drug compound was tested at two concentrations. Therefore, we have 2\*2=4 measurements corresponding to the same pair of drug compounds, each having multiple replicates for additional validation.


## Files
- `Covid-Combo.parquet`: data file (table), where each row represents morphological features obtained from Cell Painting -> CellProfiler features, per cell -> Median-averaged features, per well. First five columns contain metadata;
- `id_to_name.txt`: two-column correspondence: `batch_id` -- compound name. Note that in the data file we only use `batch_id` as a unique identifier for each measurement;
- `Covid Combo.ipynb`: python notebook with analysis based on Mahalanobis distance, i.e., ranking of combinations based on the similarity between morphological profile of a cell treatead by some drug combination and an uninfected cell;
- `mahalanobis.py`: functions and utilities for `Covid Combo.ipynb`;
- `requirements.txt`: python requirement for analysis code.


## Metadata columns
In addition to the columns corresponding to morphological features from CellProfiler, the data file  has metadata columns. Below is the description of each metadata column:
1. `batch_id`: unique identifier of each compound -- find the corresponding compound name in 'id_to_name.txt';
2. `cmpd_conc`: compound concentration in the given drug combination, in micromoles;
3. `type`: {'treated', 'infected', 'uninfected'}, i.e., drug combination, DMSO, or healthy, respectively;
4. `pair_batch_id`: batch id of another drug compound in the given combination. If not a combination (i.e., single antiviral drug, uninfected, or DMSO), then `pair_batch_id`='single';
5. `combo_id`: enumeration of unique drug combinations. If not a drug combination, then `combo_id`=-1.