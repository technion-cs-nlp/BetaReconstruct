import argparse
import os
import glob
import pandas as pd
import torch
from transformers import GPT2TokenizerFast
from transformers import ZambaForCausalLM
from transformers import AutoConfig
from ast import literal_eval
from itertools import groupby
import sys

SPECIAL_TOKENS = {
    'align': "<ALIGN>",
    'infer': "<INFER>",
    'reconstruct': "<RECONSTRUCT>"
}
OUTGROUP_SEPARATOR = "#"
SEQUENCE_SEPARATOR = "|"

def fasta_iter(fasta_name):
    """
    Given a FASTA file, yield tuples of (header, sequence).
    """
    with open(fasta_name) as fh:
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for header in faiter:
            header = next(header)[1:].strip()
            seq = "".join(s.strip() for s in next(faiter))
            yield header, seq.replace("-", "")

def parse_arguments():
    parser = argparse.ArgumentParser(description="ZAMBA inference script")
    parser.add_argument("--model_path", type=str, required=True, help="Path to the trained ZAMBA model")
    parser.add_argument("--input", type=str, required=True, help="Input CSV, FASTA file, or folder of FASTA files")
    parser.add_argument("--output", type=str, required=True, help="Output CSV, FASTA, or folder path")
    parser.add_argument("--outgroup", type=str, default=None,
                        help="Reconstruct an internal node from FASTA input: the names of the outgroup sequences, "
                             "comma-separated or as a .txt file with one name per line. All other sequences form "
                             "the ingroup, and the ancestor of the ingroup is reconstructed.")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--align", action="store_true", help="Run sequence alignment")
    group.add_argument("--infer", action="store_true", help="Run tree inference")
    group.add_argument("--reconstruct", action="store_true", help="Run ancestral sequence reconstruction")

    return parser.parse_args()

def load_model_and_tokenizer(model_path):
    print(f"Loading model from: {model_path}")
    tokenizer = GPT2TokenizerFast.from_pretrained(model_path)
    config = AutoConfig.from_pretrained(model_path)
    config.use_mamba_kernels = False

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device == "cuda":
        print("Using GPU")
        model = ZambaForCausalLM.from_pretrained(model_path)
    else:
        print("Using CPU")
        config = AutoConfig.from_pretrained(model_path)
        config.use_mamba_kernels = False
        model = ZambaForCausalLM.from_pretrained(model_path, config=config)

    model.generation_config.pad_token_id = tokenizer.pad_token_id
    model.generation_config.eos_token_id = tokenizer.eos_token_id
    model.to(device)
    model.eval()
    print("Model and tokenizer loaded successfully.")
    return model, tokenizer

def validate_species_count(seqs):
    count = len(seqs)
    print(f"Number of sequences: {count}")
    if count > 20:
        raise ValueError("Error: Too many sequences (over 20). Please reduce the number of species.")
    elif 15 <= count <= 20:
        print("Warning: Performance may degrade for >15 sequences.")

def prepare_input_string(seq_list, task_token, outgroup_list=None):
    """
    Root reconstruction, alignment and tree inference: sequences separated by '|'.
    Internal-node reconstruction (outgroup_list given): each outgroup sequence followed by '#',
    then the ingroup sequences separated by '|'.
    """
    if outgroup_list is None:
        validate_species_count(seq_list)
        return SEQUENCE_SEPARATOR.join(seq_list) + task_token
    if not seq_list:
        raise ValueError("Error: The ingroup is empty. At least one sequence must not be in the outgroup.")
    validate_species_count(list(outgroup_list) + list(seq_list))
    print(f"Internal node: {len(outgroup_list)} outgroup and {len(seq_list)} ingroup sequences")
    return "".join(seq + OUTGROUP_SEPARATOR for seq in outgroup_list) + SEQUENCE_SEPARATOR.join(seq_list) + task_token

def prepare_input_from_seq_x(seq_x, task_token):
    """
    seq_x is an internal-node input string without the task token, as in the files in data/.
    """
    parts = seq_x.split(OUTGROUP_SEPARATOR)
    outgroup_list, ingroup_list = parts[:-1], parts[-1].split(SEQUENCE_SEPARATOR)
    return prepare_input_string(ingroup_list, task_token, outgroup_list)

def read_outgroup_names(outgroup_arg):
    if os.path.isfile(outgroup_arg):
        with open(outgroup_arg) as f:
            names = [line.strip() for line in f if line.strip()]
    else:
        names = [name.strip() for name in outgroup_arg.split(",") if name.strip()]
    if not names:
        raise ValueError("Error: --outgroup is empty.")
    return names

def split_ingroup_outgroup(records, outgroup_names, path):
    """
    Split (header, sequence) records into outgroup and ingroup sequences by header.
    """
    headers = [header for header, _ in records]
    missing = [name for name in outgroup_names if name not in headers]
    if missing:
        raise ValueError(f"Error: Outgroup names not found in {path}: {', '.join(missing)}")
    outgroup_seqs = [seq for header, seq in records if header in outgroup_names]
    ingroup_seqs = [seq for header, seq in records if header not in outgroup_names]
    return ingroup_seqs, outgroup_seqs

def prepare_fasta_input(path, task, outgroup_names=None):
    records = list(fasta_iter(path))
    if outgroup_names is None:
        return prepare_input_string([seq for _, seq in records], SPECIAL_TOKENS[task])
    ingroup_seqs, outgroup_seqs = split_ingroup_outgroup(records, outgroup_names, path)
    return prepare_input_string(ingroup_seqs, SPECIAL_TOKENS[task], outgroup_seqs)

def generate_output(model, tokenizer, input_string, task):
    inputs = tokenizer(input_string, return_tensors='pt').to(model.device)
    input_len = inputs['input_ids'].shape[1]
    max_length = model.config.max_position_embeddings

    try:
        with torch.no_grad():
            outputs = model.generate(
                **inputs,
                max_length=max_length
            )
    except Exception as e:
        try:
            with torch.no_grad():
                outputs = model.generate(
                    **inputs,
                    max_new_tokens=max_length
                )
                print(f"Warning: Input too long ({input_len} tokens). Falling back to `max_new_tokens={max_length}`.")
        except Exception as e:
            raise RuntimeError(f"Model generation failed: {e}")

    decoded = tokenizer.decode(outputs[0], skip_special_tokens=False)
    decoded_tokens = "".join(decoded.split())
    parts = decoded_tokens.split(SPECIAL_TOKENS[task])
    if len(parts) >= 2:
        x, prediction = parts[0], parts[1]
    else:
        raise RuntimeError(f"Failed to process output: {decoded_tokens}")
    if tokenizer.decode(tokenizer.eos_token_id) in prediction:
        prediction = prediction.split(tokenizer.decode(tokenizer.eos_token_id))[0]
    if "!" in prediction:
        prediction = prediction.replace("!", ";")
    if "*" in prediction:
        prediction = prediction.replace("*", "-")
    return prediction

def prepare_csv_row_input(row, columns, task):
    """
    Internal-node reconstruction uses the seq_x column, or else the outgroup and ingroup columns.
    Otherwise all sequences in unaligned_seqs are used (root reconstruction, alignment, tree inference).
    """
    if task == 'reconstruct' and 'seq_x' in columns:
        return prepare_input_from_seq_x(row['seq_x'], SPECIAL_TOKENS[task])
    if task == 'reconstruct' and 'outgroup' in columns and 'ingroup' in columns:
        return prepare_input_string(literal_eval(row['ingroup']), SPECIAL_TOKENS[task], literal_eval(row['outgroup']))
    seqs = literal_eval(row['unaligned_seqs'])
    return prepare_input_string(seqs, SPECIAL_TOKENS[task])

def handle_csv_input(model, tokenizer, path, output_path, task):
    print(f"Processing CSV: {path}")
    df = pd.read_csv(path)
    internal_node_columns = 'seq_x' in df.columns or ('outgroup' in df.columns and 'ingroup' in df.columns)
    if 'unaligned_seqs' not in df.columns and not (task == 'reconstruct' and internal_node_columns):
        raise ValueError("CSV must contain a column named 'unaligned_seqs' (or, to reconstruct internal nodes, "
                         "'seq_x' or 'outgroup' and 'ingroup').")
    if task == 'reconstruct':
        if internal_node_columns:
            print("Reconstructing internal nodes (the ancestor of the ingroup).")
        else:
            print("No outgroup columns found: reconstructing the root of all sequences.")
    if f'prediction_{task}' not in df.columns:
        df[f'prediction_{task}'] = None

    for i, row in df.iterrows():
        if isinstance(row[f'prediction_{task}'], str):
            continue
        try:
            prompt = prepare_csv_row_input(row, df.columns, task)
            result = generate_output(model, tokenizer, prompt, task)
            print(f"Result: {result}")
            df.at[i, f'prediction_{task}'] = result
            print(f"Processed row {i+1}/{len(df)}")
        except Exception as e:
            print(f"Error processing row {i}: {e}")

    df.to_csv(output_path, index=False)
    print(f"Saved predictions to: {output_path}")

def handle_fasta_input(model, tokenizer, path, output_path, task, outgroup_names=None):
    print(f"Processing FASTA: {path}")
    prompt = prepare_fasta_input(path, task, outgroup_names)
    result = generate_output(model, tokenizer, prompt, task)
    print(f"Result: {result}")
    with open(output_path, 'w') as f:
        f.write(">prediction\n")
        f.write(result + "\n")
    print(f"Saved prediction to: {output_path}")

def find_ingroup_outgroup_folders(folder_path):
    """
    Subfolders that contain ingroup.fasta and outgroup.fasta, as in data/mammalian_validation.
    """
    return sorted(
        sub for sub in glob.glob(os.path.join(folder_path, "*"))
        if os.path.isfile(os.path.join(sub, "ingroup.fasta")) and os.path.isfile(os.path.join(sub, "outgroup.fasta"))
    )

def handle_fasta_folder(model, tokenizer, folder_path, output_folder, task, outgroup_names=None):
    print(f"Processing folder of FASTA files: {folder_path}")
    os.makedirs(output_folder, exist_ok=True)
    fasta_files = sorted(glob.glob(os.path.join(folder_path, "*.fasta")))
    group_folders = find_ingroup_outgroup_folders(folder_path) if task == 'reconstruct' else []
    if not fasta_files and not group_folders:
        raise FileNotFoundError("No .fasta files (or subfolders with ingroup.fasta and outgroup.fasta) found in the specified folder.")
    for fasta_path in fasta_files:
        try:
            prompt = prepare_fasta_input(fasta_path, task, outgroup_names)
            result = generate_output(model, tokenizer, prompt, task)
            print(f"Result: {result}")
            out_path = os.path.join(output_folder, os.path.basename(fasta_path))
            with open(out_path, 'w') as f:
                f.write(">prediction\n")
                f.write(result + "\n")
            print(f"Saved: {out_path}")
        except Exception as e:
            print(f"Failed on {fasta_path}: {e}")
    for group_folder in group_folders:
        try:
            ingroup_seqs = [seq for _, seq in fasta_iter(os.path.join(group_folder, "ingroup.fasta"))]
            outgroup_seqs = [seq for _, seq in fasta_iter(os.path.join(group_folder, "outgroup.fasta"))]
            prompt = prepare_input_string(ingroup_seqs, SPECIAL_TOKENS[task], outgroup_seqs)
            result = generate_output(model, tokenizer, prompt, task)
            print(f"Result: {result}")
            out_path = os.path.join(output_folder, os.path.basename(os.path.normpath(group_folder)) + ".fasta")
            with open(out_path, 'w') as f:
                f.write(">prediction\n")
                f.write(result + "\n")
            print(f"Saved: {out_path}")
        except Exception as e:
            print(f"Failed on {group_folder}: {e}")

def main(args):
    print(args)
    if args.align:
        task = 'align'
    elif args.infer:
        task = 'infer'
    elif args.reconstruct:
        task = 'reconstruct'
    else:
        raise ValueError("Error: You must specify one task: --align, --infer, or --reconstruct")

    outgroup_names = None
    if args.outgroup is not None:
        if task != 'reconstruct':
            raise ValueError("Error: --outgroup can only be used with --reconstruct.")
        outgroup_names = read_outgroup_names(args.outgroup)

    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input path does not exist: {args.input}")

    model, tokenizer = load_model_and_tokenizer(args.model_path)

    if os.path.isdir(args.input):
        handle_fasta_folder(model, tokenizer, args.input, args.output, task, outgroup_names)
    elif args.input.endswith(".csv"):
        if outgroup_names is not None:
            print("Warning: --outgroup is ignored for CSV input; use the seq_x column or the outgroup and ingroup columns.")
        handle_csv_input(model, tokenizer, args.input, args.output, task)
    elif args.input.endswith((".fasta", ".fa")):
        handle_fasta_input(model, tokenizer, args.input, args.output, task, outgroup_names)
    else:
        raise ValueError("Unsupported input type. Use a .csv, .fasta file, or folder of .fasta files.")

if __name__ == "__main__":
    main(parse_arguments())
