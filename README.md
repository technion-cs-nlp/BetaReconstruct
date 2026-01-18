# Ancestral sequence reconstruction using generative models

## Abstract:
Ancestral sequence reconstruction (ASR) is a foundational task in evolutionary biology, providing insights into the molecular past and guiding studies of protein function and adaptation. Conventional ASR methods rely on a multiple sequence alignment (MSA), a phylogenetic tree, and an evolutionary model. However, the underlying alignments and trees are often uncertain, and existing models typically focus on substitutions and do not explicitly account for insertion-deletion (indel) processes. Here, we introduce BetaReconstruct, a novel generative approach to ASR that harnesses recent advances in natural language processing (NLP) and hybrid transformer architectures. Our model was initially trained on large-scale simulated datasets with gold-standard ancestral sequences and subsequently on real-world protein sequences. The reconstruction requires neither MSAs nor phylogenetic trees. We demonstrate that BetaReconstruct generalizes robustly across diverse evolutionary scenarios and reconstructs ancestral sequences more accurately than maximum-likelihood-based pipelines. We additionally provide evidence that the generative-model ASR approach is also more accurate when analyzing empirical datasets. This work provides a scalable, alignment-free strategy for ASR and highlights the ability of data-driven models to capture evolutionary signals beyond the reach of traditional methods.

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
 - [Research paper](https://academic.oup.com/bioinformatics/article/40/4/btae196/7645044): For detailed methodology and performance benchmarks, please see our paper: *"Effect of tokenization on transformers for biological sequences"* (Dotan et al., 2024; Bioinformatics).


## Align-Infer-Reconstruct pipeline

This script runs a ZAMBA-based language model for performing one of three tasks in molecular phylogenetics:
 - Alignment of unaligned sequences
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

export MODEL_PATH="dotan1111/BetaReconstruct_Configuratio2"
export INPUT="$BETA_RECONSTRUCT/example_inputs_for_align_infer_or_reconstruct/example1.fasta"
export OUTPUT="$BETA_RECONSTRUCT/outputs/example1.fasta"

cd $BETA_RECONSTRUCT

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
export OUTPUT="$BETA_RECONSTRUCT/outputs/results_folder"

cd $BETA_RECONSTRUCT

python "align_infer_or_reconstruct.py" \
    --model_path $MODEL_PATH \
    --input $INPUT \
    --output $OUTPUT \
    --align
```



## Inferring ancestral sequence for many species
This Python script performs hierarchical ancestral sequence reconstruction using a language model (ZambaForCausalLM) in a multi-level strategy. Sequences are grouped based on clade and progressively merged and reconstructed level by level until a root ancestral sequence is generated. We note that the models were fine-tuned on **Mammals data** only (OrthoMaM; Allio et al., 2024; Nucleic Acids Research).

### Overview
The method uses a grouping strategy to divide sequences into manageable groups, generating intermediate ancestral sequences at each level using a causal language model. This continues hierarchically until a single root ancestor is reconstructed.

The script supports FASTA files as input and expects species groupings in text format. Output is saved in JSON format containing group information and the final prediction.

#### Input Requirements
1. FASTA Files
 - Each fasta file should contain multiple sequences for one gene or family.
 - Filenames should follow the format FAMILY_NAME.fasta.

2. Species List Folder
 - Contains .txt files, each listing species names (one per line) for a clade (e.g., an order or family).
- Used to map sequences to their respective taxonomic groups.

#### Input Arguments
Argument	Description
+ --model_path	(str, Required) Path or name of the pretrained ZambaForCausalLM model.
+ --fasta_files_folder	(str, Required) Directory containing input *.fasta* files. Filenames must match species in species_list_folder.
+ --results_folder	(str, Required) Directory where json result files will be saved.
+ --species_list_folder	(str, Required) Folder containing .txt files of species grouped by clades (e.g., orders). Used to divide input data.
+ --model_name	(str, Optional) Tag added to result filenames (e.g., for model versioning or experiment tracking).

#### Output
For each *fasta file*, a corresponding *json* is saved in the results_folder, containing:

Reconstruction groups at each hierarchical level

The final root ancestral sequence prediction

#### Examples

```
export BETA_RECONSTRUCT="<PATH_TO_BETA_RECONSTRUCT>"
conda activate $BETA_RECONSTRUCT/python_env/

export HF_DATASETS_CACHE="$BETA_RECONSTRUCT/python_env/cache/"
export HF_HOME="$BETA_RECONSTRUCT/python_env/cache/"

export MODEL_PATH="dotan1111/BetaReconstruct_Mammals_Configuration1"
export FASTA_FILES_FOLDER="$BETA_RECONSTRUCT/fasta_val_files/"
export RESULTS_FOLDER="$BETA_RECONSTRUCT/results/"
export SPEICEIS_LIST_FOLDER="$BETA_RECONSTRUCT/species_lists/"

cd $BETA_RECONSTRUCT

python "predict_ancestral_hierarchical_approach.py" \
    --model_path $MODEL_PATH \
    --fasta_files_folder $FASTA_FILES_FOLDER \
    --results_folder $RESULTS_FOLDER \
    --species_list_folder $SPEICEIS_LIST_FOLDER

```

