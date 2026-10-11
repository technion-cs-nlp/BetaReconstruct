# Simulated test datasets

This folder contains the simulated test sets used to evaluate BetaReconstruct in *Ancestral sequence reconstruction using generative models* (Dotan et al.). They are provided so that other ASR methods can be evaluated on exactly the same instances.

In every instance, the task is to reconstruct the ancestral sequence at an internal node of the tree (the target node). The input is a set of unaligned sequences: the ingroup sequences, which are the descendants of the target node, and two outgroup sequences. The 'true' ancestral sequence at the target node is known from the simulation.

## Files

| File | Dataset in the paper | Instances | Sequences per instance | Results in the paper |
|---|---|---|---|---|
| `DI3_internal_nodes_test.csv` | DI3 | 1,000 (200 per species count) | 10–14 (8–12 ingroup + 2 outgroup) | Table 1, Table S2 |
| `DMP_multi_partition_test.csv` | DMP (spatially heterogeneous model) | 1,000 (200 per species count) | 10–14 (8–12 ingroup + 2 outgroup) | Table 2 |
| `OOD_replacement_matrix_test.csv` | Out-of-distribution: replacement matrices | 800 (200 each: WAG, Dayhoff, JTT, LG) | 14 (12 ingroup + 2 outgroup) | Fig. 2A, Table S5a |
| `OOD_indel_length_test.csv` | Out-of-distribution: indel length distributions | 600 (200 each: Zipfian, geometric, Poisson) | 14 (12 ingroup + 2 outgroup) | Fig. 2B, Table S5b |

In the out-of-distribution files, WAG and Zipfian are the in-distribution settings (BetaReconstruct was trained on WAG and Zipfian data); the other settings were not seen during training. The simulation parameters of each dataset are given in Supplementary Information S1 and Table S1 of the paper.

## Columns

| Column | Description |
|---|---|
| `instance_id` | Index of the instance within the file (0-based). Use it to match predictions to instances. |
| `num_species` | Total number of sequences in the instance (ingroup + 2 outgroup). |
| `unaligned_seqs` | All unaligned sequences, as a Python list: the two outgroup sequences first, then the ingroup sequences. |
| `outgroup` | The two outgroup sequences (Python list). |
| `ingroup` | The ingroup sequences, i.e., the descendants of the target node (Python list). |
| `seq_x` | The input string given to BetaReconstruct: each outgroup sequence followed by `#`, then the ingroup sequences separated by `\|`. The model receives `seq_x` followed by the `<RECONSTRUCT>` token. |
| `true` | The 'true' (simulated) ancestral sequence at the target node, unaligned. |
| `tree_true` | `DI3_internal_nodes_test.csv` only. The true tree in Newick format. Leaves are named `l0` … `l{N-1}`, where N = `num_species`. Ingroup sequence k is leaf `l{k}`, and the two outgroup sequences are leaves `l{N-2}` and `l{N-1}`. The target node is the child of the root that contains the ingroup leaves. |
| `sub_model` | `OOD_replacement_matrix_test.csv` only. Replacement matrix used to simulate the instance: `WAG`, `DAYHOFF`, `JTT`, or `LG`. |
| `indel_model` | `OOD_indel_length_test.csv` only. Indel length distribution used to simulate the instance: `ZIP` (Zipfian), `GEO` (geometric), or `POI` (Poisson). |

## Loading the data

```python
import pandas as pd
from ast import literal_eval

df = pd.read_csv("DI3_internal_nodes_test.csv")
for col in ["unaligned_seqs", "outgroup", "ingroup"]:
    df[col] = df[col].apply(literal_eval)
```

## Evaluation

In the paper, each reconstructed sequence was globally aligned to the `true` ancestral sequence, and the fractions of matches, mismatches, false ancestral insertions, and false ancestral deletions were computed, together with the Levenshtein distance (the sum of the three error types). All metrics are normalized by the length of the true ancestral sequence (see Methods, "Evaluating performance").
