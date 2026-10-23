# Ancestral sequence reconstruction using generative models

## Abstract:
Ancestral sequence reconstruction (ASR) is a foundational task in evolutionary biology, providing insights into the molecular past and guiding studies of protein function and adaptation. Conventional ASR methods rely on a multiple sequence alignment (MSA), a phylogenetic tree, and an evolutionary model. However, the underlying alignments and trees are often uncertain, and existing models typically focus on substitutions and do not explicitly account for insertion-deletion (indel) processes. Here, we introduce BetaReconstruct, a novel generative approach to ASR that harnesses recent advances in natural language processing (NLP) and hybrid transformer architectures. Our model was initially trained on large-scale simulated datasets with gold-standard ancestral sequences and subsequently on real-world protein sequences. The reconstruction requires neither MSAs nor phylogenetic trees; to reconstruct an internal node, the user only specifies which sequences descend from it (ingroup) and which do not (outgroup). BetaReconstruct achieved accuracy comparable to that of maximum-likelihood and Bayesian methods in reconstructing ancestral sequences at internal nodes, and inferred fewer false ancestral insertions. Its accuracy remained similar under evolutionary models not seen during training. When the evolutionary process varied along the sequence, BetaReconstruct was in some cases more accurate than maximum-likelihood methods. . This work provides a scalable, alignment-free strategy for ASR and highlights the ability of data-driven models to capture evolutionary signals beyond the reach of traditional methods.

![alt text](https://github.com/technion-cs-nlp/BetaReconstruct/blob/main/outline_image.png)
Illustration of ASR prediction using BetaReconstruct. (a): Consider the evolution of the sequence “AAMM”, to the proteins: “AAM”, “AYM” and “ATMMM”. (b) BetaReconstruct pipeline. (Ⅰ): The input to the model, unaligned protein sequences; (Ⅱ): The protein sequences are concatenated with special characters between; (Ⅲ): The model processes the input and predicts the output, which is the ancestral sequence (Ⅳ).

## Public models 
We have released six pre-trained models on the HuggingFace Hub, available for public download. These models are categorized into two main groups based on their training objectives and datasets.
1. **Mammalian Ancestral Reconstruction Models**

These models were specifically fine-tuned on mammalian protein sequences to optimize accuracy for Ancestral Sequence Reconstruction within mammalian lineages.
- [BetaReconstruct_Mammals_Configuration1](https://huggingface.co/dotan1111/BetaReconstruct_Mammals_Configuration1)
- [BetaReconstruct_Mammals_Configuration2](https://huggingface.co/dotan1111/BetaReconstruct_Mammals_Configuration2)
- [BetaReconstruct_Mammals_Configuration3](https://huggingface.co/dotan1111/BetaReconstruct_Mammals_Configuration3)

2. **General Simulated Data Models**

These models were trained on large-scale simulated datasets. They are multi-purpose tools designed to generate:
+ Multiple Sequence Alignments
+ Phylogenetic Trees
+ Ancestral Sequences

Links to models:
 - [BetaReconstruct_Configuration1](https://huggingface.co/dotan1111/BetaReconstruct_Configuration1)
 - [BetaReconstruct_Configuration2](https://huggingface.co/dotan1111/BetaReconstruct_Configuration2)
 - [BetaReconstruct_Configuration3](https://huggingface.co/dotan1111/BetaReconstruct_Configuration3)

**[IMPORTANT] Note on Generalization**: While we have verified the models' ability to generalize to out-of-distribution data (as detailed in the main text), performance is highest when inference is performed on data distributions similar to those used during training.


### Tokenization
The models utilize a data-driven tokenizer specifically optimized for evolutionary related protein sequences. Instead of relying on standard single-amino-acid characters, this approach allows the model to learn and represent frequently occurring sequence patterns as distinct tokens. For a comprehensive technical breakdown of the tokenizer architecture and its training process, please refer to the main text of our publication.

Related Resources
 - [BiologicalTokenizers repository](https://github.com/technion-cs-nlp/BiologicalTokenizers): You can find the source code and training scripts for biological tokenizers.
 - [Research paper](https://academic.oup.com/bioinformatics/article/40/4/btae196/7645044): For detailed methodology and performance benchmarks, please see our paper: *"Effect of tokenization on transformers for biological sequences"* (2024; Bioinformatics).


## Align-Infer-Reconstruct pipeline

This script runs a ZAMBA-based language model for performing one of three tasks in molecular phylogenetics:
 - Alignment of unaligned sequences: See also: [Multiple sequence alignment as a sequence-to-sequence learning problem (2023; ICLR)](https://openreview.net/forum?id=8efJYMBrNb) and [BetaAlign: a deep learning approach for multiple sequence alignment (2025; Bioinformatics)](https://academic.oup.com/bioinformatics/article/41/1/btaf009/7945664)
 - Inference of phylogenetic relationships
 - Reconstruction of ancestral sequences

The model works with FASTA or CSV inputs.
The models were trained on simulated data, containing 10 - 14 species.
Input examples are provided in example_inputs_for_align_infer_or_reconstruct

### Input Requirements
Argument	Description
+ --model_path (str, Required) Path or name of the pretrained ZambaForCausalLM model.
+ --align (flag) Use this flag to perform sequence alignment.
+ --infer (flag) Use this flag to perform phylogenetic inference.
+ --reconstruct (flag) Use this flag to reconstruct ancestral sequences.
+ --input (str, Required) Path to input file or folder. Must be either:
  + a. a csv file with unaligned_seqs and num_species columns
  + b. a fasta file
  + c. a folder containing .fasta files
+ --output (str, Required) Path to the output file (for CSV or FASTA) or folder (for batch FASTA processing).

#### Examples

##### Reconstruct an internal node using a fasta file and an outgroup
The sequences named in `--outgroup` (a .txt file with one name per line, or comma-separated names) form the outgroup, and all other sequences form the ingroup. The ancestor of the ingroup is reconstructed.
```
export BETA_RECONSTRUCT="<PATH_TO_BETA_RECONSTRUCT>"
conda activate $BETA_RECONSTRUCT/python_env/
export HF_DATASETS_CACHE="$BETA_RECONSTRUCT/python_env/cache/"
export HF_HOME="$BETA_RECONSTRUCT/python_env/cache/"
export MODEL_PATH="<INTERNAL_NODE_MODEL>"
export INPUT="$BETA_RECONSTRUCT/example_inputs_for_align_infer_or_reconstruct/internal_node_example.fasta"
export OUTGROUP="$BETA_RECONSTRUCT/example_inputs_for_align_infer_or_reconstruct/internal_node_outgroup.txt"
export OUTPUT="$BETA_RECONSTRUCT/outputs/internal_node_example_prediction.fasta"
cd $BETA_RECONSTRUCT
mkdir -p "$(dirname "$OUTPUT")"
python "align_infer_or_reconstruct.py" \
    --model_path $MODEL_PATH \
    --input $INPUT \
    --output $OUTPUT \
    --outgroup $OUTGROUP \
    --reconstruct
```
##### Reconstruct internal nodes using a csv file
Each row is one instance. The script uses the `seq_x` column (or the `outgroup` and `ingroup` columns) and adds a `prediction_reconstruct` column. The test sets in `data/` can be used the same way.
```
export INPUT="$BETA_RECONSTRUCT/example_inputs_for_align_infer_or_reconstruct/internal_node_test.csv"
export OUTPUT="$BETA_RECONSTRUCT/outputs/internal_node_predictions.csv"
python "align_infer_or_reconstruct.py" \
    --model_path $MODEL_PATH \
    --input $INPUT \
    --output $OUTPUT \
    --reconstruct
```

##### Phylogenetic tree inference using a csv file
```
export BETA_RECONSTRUCT="<PATH_TO_BETA_RECONSTRUCT>"
conda activate $BETA_RECONSTRUCT/python_env/

export HF_DATASETS_CACHE="$BETA_RECONSTRUCT/python_env/cache/"
export HF_HOME="$BETA_RECONSTRUCT/python_env/cache/"

export MODEL_PATH="dotan1111/BetaReconstruct_Configuration1"
export INPUT="$BETA_RECONSTRUCT/example_inputs_for_align_infer_or_reconstruct/test.csv"
export OUTPUT="$BETA_RECONSTRUCT/outputs/trees.csv"

cd $BETA_RECONSTRUCT
mkdir -p "$(dirname "$OUTPUT")"

python "align_infer_or_reconstruct.py" \
    --model_path $MODEL_PATH \
    --input $INPUT \
    --output $OUTPUT \
    --infer
```

##### Align using a fasta file
```
export BETA_RECONSTRUCT="<PATH_TO_BETA_RECONSTRUCT>"
conda activate $BETA_RECONSTRUCT/python_env/

export HF_DATASETS_CACHE="$BETA_RECONSTRUCT/python_env/cache/"
export HF_HOME="$BETA_RECONSTRUCT/python_env/cache/"

export MODEL_PATH="dotan1111/BetaReconstruct_Configuration2"
export INPUT="$BETA_RECONSTRUCT/example_inputs_for_align_infer_or_reconstruct/example1.fasta"
export OUTPUT="$BETA_RECONSTRUCT/outputs/example1_aligned.fasta"

cd $BETA_RECONSTRUCT
mkdir -p "$(dirname "$OUTPUT")"

python "align_infer_or_reconstruct.py" \
    --model_path $MODEL_PATH \
    --input $INPUT \
    --output $OUTPUT \
    --align
```

##### Reconstruct using a folder
```
export BETA_RECONSTRUCT="<PATH_TO_BETA_RECONSTRUCT>"
conda activate $BETA_RECONSTRUCT/python_env/

export HF_DATASETS_CACHE="$BETA_RECONSTRUCT/python_env/cache/"
export HF_HOME="$BETA_RECONSTRUCT/python_env/cache/"

export MODEL_PATH="dotan1111/BetaReconstruct_Configuration3"
export INPUT="$BETA_RECONSTRUCT/example_inputs_for_align_infer_or_reconstruct/folder_with_fasta"
export OUTPUT="$BETA_RECONSTRUCT/outputs/ASR_results_folder"

cd $BETA_RECONSTRUCT
mkdir -p "$OUTPUT"

python "align_infer_or_reconstruct.py" \
    --model_path $MODEL_PATH \
    --input $INPUT \
    --output $OUTPUT \
    --reconstruct
```



### Mammalian proteins (OrthoMaM)
The mammalian test set is in [`data/mammalian_validation/`](data/mammalian_validation/): 40 protein families from OrthoMaM version 12 (Allio et al., 2024), the mammalian test families of the paper. In each family, the task is to reconstruct the ancestor of the placental mammals:
 - **Ingroup:** the placental species (Placentalia), about 164 per family.
 - **Outgroup:** the marsupial species (orders Dasyuromorphia, Didelphimorphia, Diprotodontia, Microbiotheria (*Dromiciops gliroides*) and Peramelemorphia), up to 8 per family.
 - The monotremes (platypus and echidna), which are neither placental nor marsupial, were removed.
Each family has its own folder with `ingroup.fasta` and `outgroup.fasta` (unaligned sequences, named by species). The folder also contains `mammalian_test_families.csv` (all families in the input format of the simulated test sets), `family_counts.csv` (number of ingroup and outgroup species per family), and the script that created these files from OrthoMaM. Families 283297 and 5535 contain no marsupials, so their ingroup ancestor cannot be reconstructed. In the paper, the reconstructed ancestors were evaluated through their predicted structures (Supplementary Information S5).
Because each family has about 170 species, the ancestors are reconstructed hierarchically with `predict_ancestral_hierarchical_approach.py` (see [Inferring ancestral sequence for many species](#inferring-ancestral-sequence-for-many-species)), which reads the folder directly.

#### Running example
Reconstruct the ancestor of the placental mammals for one family (2331), then for all families:
```
export BETA_RECONSTRUCT="<PATH_TO_BETA_RECONSTRUCT>"
conda activate $BETA_RECONSTRUCT/python_env/
export HF_DATASETS_CACHE="$BETA_RECONSTRUCT/python_env/cache/"
export HF_HOME="$BETA_RECONSTRUCT/python_env/cache/"
export MODEL_PATH="<MAMMALIAN_INTERNAL_NODE_MODEL>"
export RESULTS_FOLDER="$BETA_RECONSTRUCT/results/mammalian/"
export SPECIES_LIST_FOLDER="$BETA_RECONSTRUCT/species_lists/"
cd $BETA_RECONSTRUCT
mkdir -p "$RESULTS_FOLDER"
# One family
python "predict_ancestral_hierarchical_approach.py" \
    --model_path $MODEL_PATH \
    --fasta_files_folder "$BETA_RECONSTRUCT/data/mammalian_validation/2331_NT_AL_AA/" \
    --results_folder $RESULTS_FOLDER \
    --species_list_folder $SPECIES_LIST_FOLDER
# All families (283297 and 5535 are skipped: no outgroup)
python "predict_ancestral_hierarchical_approach.py" \
    --model_path $MODEL_PATH \
    --fasta_files_folder "$BETA_RECONSTRUCT/data/mammalian_validation/" \
    --results_folder $RESULTS_FOLDER \
    --species_list_folder $SPECIES_LIST_FOLDER
```
Each family gives `results/mammalian/<family>.json`. The reconstructed placental ancestor is in its `prediction` field:
```
python -c "import json; print(json.load(open('results/mammalian/2331_NT_AL_AA.json'))['prediction'])"
```

## APA

```
Edo Dotan, Elya Wygoda, Asaf Schers, Iris Lyubman, Yonatan Belinkov and Tal Pupko, Ancestral sequence reconstruction using generative model, bioRxiv 2026; https://www.biorxiv.org/content/early/2026/01/21/2026.01.18.700141
```


## BibTeX
```
@article{dotan_ancestral_2026,
	title = {Ancestral sequence reconstruction using generative models},
	url = {https://www.biorxiv.org/content/early/2026/01/21/2026.01.18.700141},
	doi = {10.64898/2026.01.18.700141},
	abstract = {Ancestral sequence reconstruction (ASR) is a foundational task in evolutionary biology, providing insights into the molecular past and guiding studies of protein function and adaptation. Conventional ASR methods rely on a multiple sequence alignment (MSA), a phylogenetic tree, and an evolutionary model. However, the underlying alignments and trees are often uncertain, and existing models typically focus on substitutions and do not explicitly account for insertion-deletion (indel) processes. Here, we introduce BetaReconstruct, a novel generative approach to ASR that harnesses recent advances in natural language processing (NLP) and hybrid transformer architectures. Our model was initially trained on large-scale simulated datasets with gold-standard ancestral sequences and subsequently on real-world protein sequences. The reconstruction requires neither MSAs nor phylogenetic trees. We demonstrate that BetaReconstruct generalizes robustly across diverse evolutionary scenarios and reconstructs ancestral sequences more accurately than maximum-likelihood-based pipelines. We additionally provide evidence that the generative-model ASR approach is also more accurate when analyzing empirical datasets. This work provides a scalable, alignment-free strategy for ASR and highlights the ability of data-driven models to capture evolutionary signals beyond the reach of traditional methods.},
	journal = {bioRxiv},
	author = {Dotan, Edo and Wygoda, Elya and Schers, Asaf and Lyubman, Iris and Belinkov, Yonatan and Pupko, Tal},
	year = {2026},
	note = {Publisher: Cold Spring Harbor Laboratory
\_eprint: https://www.biorxiv.org/content/early/2026/01/21/2026.01.18.700141.full.pdf},
}
```
