# Reconstructed ancestral sequences

This folder contains the ancestral sequences reconstructed by each method on the test sets of *Ancestral sequence reconstruction using generative models* (Dotan et al.). They are the reconstructions behind the results reported in the paper, provided so that the results can be checked and other ASR methods can be compared on exactly the same instances.

The methods are BetaReconstruct, FastML (Pupko et al., 2000), ARPIP (Jowkar et al., 2023) and BAli-Phy (Redelings, 2021). All methods received the same unaligned sequences. FastML and ARPIP were run on a MAFFT alignment and a maximum-likelihood tree (IQ-TREE) rooted with the outgroup, under the WAG model; see Methods in the paper for the full settings.

## Files

Each file has the same rows and input columns as the matching test set in [`data/`](../data/) (`<dataset>_test.csv` → `<dataset>_preds.csv`), with the predictions added as columns. The mammalian file has one row per family in [`data/mammalian_validation/`](../data/mammalian_validation/).

| File | Dataset in the paper | Rows | BetaReconstruct | FastML | ARPIP | BAli-Phy | Results in the paper |
|---|---|---|---|---|---|---|---|
| `DI3_internal_nodes_preds.csv` | DI3 | 1,000 | 1,000 | 987 | 987 | 999 | Table 1 |
| `DMP_multi_partition_preds.csv` | DMP (spatially heterogeneous model) | 1,000 | 1,000 | 983 | 983 | 1,000 | Table 2 |
| `OOD_replacement_matrix_preds.csv` | Out-of-distribution: replacement matrices | 800 | 800 | 789 | 789 | 798 | Fig. 2A |
| `OOD_indel_length_preds.csv` | Out-of-distribution: indel length distributions | 600 | 600 | 589 | 589 | 599 | Fig. 2B |
| `mammalian_test_families_preds.csv` | OrthoMaM mammalian families | 40 | 37 | 30 | 30 | – | Table 3 |

The method columns give the number of rows with a reconstruction. An empty cell means that the method produced no reconstruction for that instance. BAli-Phy was not run on the mammalian families, because it did not converge on families of this size.

## Columns

### Simulated test sets

| Column | Description |
|---|---|
| `instance_id` | Index of the instance; identical to `instance_id` in the matching file in `data/`. |
| `num_species` | Total number of sequences in the instance (ingroup + 2 outgroup). |
| `outgroup` | The two outgroup sequences (Python list). |
| `ingroup` | The ingroup sequences, i.e., the descendants of the target node (Python list). |
| `seq_x` | The input string given to BetaReconstruct: each outgroup sequence followed by `#`, then the ingroup sequences separated by `\|`. |
| `true` | The 'true' (simulated) ancestral sequence at the target node, unaligned. |
| `tree_true` | `DI3_internal_nodes_preds.csv` only. The true tree in Newick format (see `data/README.md`). |
| `sub_model` | `OOD_replacement_matrix_preds.csv` only. Replacement matrix used to simulate the instance: `WAG`, `DAYHOFF`, `JTT`, or `LG`. |
| `indel_model` | `OOD_indel_length_preds.csv` only. Indel length distribution used to simulate the instance: `ZIP` (Zipfian), `GEO` (geometric), or `POI` (Poisson). |
| `pred_betareconstruct` | Ancestral sequence reconstructed by BetaReconstruct. |
| `pred_fastml` | Ancestral sequence reconstructed by FastML. |
| `pred_arpip` | Ancestral sequence reconstructed by ARPIP. |
| `pred_baliphy` | Ancestral sequence reconstructed by BAli-Phy. |

All reconstructed sequences are unaligned, and are the reconstruction at the target node: the ancestor of the ingroup sequences.

### Mammalian families

`mammalian_test_families_preds.csv` has one row per OrthoMaM family. The ingroup is the placental species and the outgroup the marsupial species (see the "Mammalian proteins (OrthoMaM)" section of the main README).

| Column | Description |
|---|---|
| `family` | OrthoMaM file name (NCBI Gene ID + `_NT_AL_AA`); the matching folder in `data/mammalian_validation/`. |
| `validation_file` | The OrthoMaM FASTA file the row was built from. |
| `outgroup`, `ingroup` | The outgroup and ingroup sequences, unaligned (Python lists). |
| `outgroup_names`, `ingroup_names` | Species names, in the order of `outgroup` and `ingroup`. |
| `num_species` | Number of ingroup and outgroup sequences. |
| `seq_x` | The sequences in the input format of BetaReconstruct (each outgroup sequence followed by `#`, then the ingroup sequences separated by `\|`). Because the families are large, BetaReconstruct reconstructed them hierarchically (`predict_ancestral_hierarchical_approach.py`), not from this single string. |
| `true`, `tree_true` | Empty: the true ancestor and tree are unknown for real proteins. |
| `pred_betareconstruct`, `pred_fastml`, `pred_arpip` | The ancestor of the ingroup (the placental mammals) reconstructed by each method. |

## Evaluation

In the paper, each reconstructed sequence of the simulated test sets was globally aligned to the `true` ancestral sequence, and the fractions of matches, mismatches, false ancestral insertions, and false ancestral deletions were computed, together with the Levenshtein distance (the sum of the three error types). All metrics are normalized by the length of the true ancestral sequence (see Methods, "Evaluating performance"). The mammalian reconstructions were evaluated through their predicted structures (Supplementary Information S5).

## Loading the data

```python
import pandas as pd
from ast import literal_eval

df = pd.read_csv("DI3_internal_nodes_preds.csv")
for col in ["outgroup", "ingroup"]:
    df[col] = df[col].apply(literal_eval)
methods = ["pred_betareconstruct", "pred_fastml", "pred_arpip", "pred_baliphy"]
complete = df.dropna(subset=methods)  # instances reconstructed by all methods
```

