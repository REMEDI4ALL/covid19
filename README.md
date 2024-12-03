# Covid: Drug Combinations

This folder contains Cell Painting results for drug combinations against COVID-19, including some analysis code.

Data consist of the CellProfiler features that were first obtained for each cell and subsequently averaged by median within each well. Among the resulting averaged features, we have three main types of measurements (cells):
1. `treated`: features of cells that were treated either with a two-drug combination or a single antiviral drug for control (Nirmatrelvir, GC376, Remdesivir);
2. `infected`: features of cells that were exposed to Covid and no compound except the solvent (DMSO) was added;
3. `uninfected`: features of healthy cells.

## Experiment
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
- `data/Covid-Combo.parquet.gzip`: data file (table), where each row represents morphological features obtained from Cell Painting -> CellProfiler features, per cell -> Median-averaged features, per well. First five columns contain metadata;
- `data/id_to_name.txt`: two-column correspondence {batch_id: compound name}. Note that in the data file we only use batch_id as a unique identifier for each measurement;
- `Covid Combo.ipynb`: python notebook with analysis based on Mahalanobis distance, i.e., ranking of combinations based on the similarity between morphological profile of a cell treatead by some drug combination and uninfected cells;
- `mahalanobis.py`: functions and utilities for `Covid Combo.ipynb`;
- `requirements.txt`: python requirement for analysis code.


## Metadata columns
In addition to the columns corresponding to morphological features from CellProfiler, we have metadata columns:
1. `type`: {'treated', 'infected', 'uninfected'}, i.e., drug (either combination or single control), DMSO, and healthy, respectively;
2. `batch_id`: unique identifier of the first drug in the given combination (find the corresponding compound name in 'id_to_name.txt');
3. `pair_batch_id`: id of the second drug in the given combination. If not a combination, i.e., single antiviral control, uninfected, or DMSO, then `pair_batch_id`='single';
4. `conc`: concentration of the first drug in the given combination, in micromoles;
5. `pair_conc`: concentration of the second drug in the given combination, in micromoles. If not a combination, i.e., single antiviral control, uninfected, or DMSO, then `pair_conc`=-1;