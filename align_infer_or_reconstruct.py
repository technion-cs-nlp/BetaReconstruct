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

def prepare_input_string(seq_list, task_token):
    validate_species_count(seq_list)
    return "|".join(seq_list) + task_token

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

def handle_csv_input(model, tokenizer, path, output_path, task):
    print(f"Processing CSV: {path}")
    df = pd.read_csv(path)
    if 'unaligned_seqs' not in df.columns:
        raise ValueError("CSV must contain a column named 'unaligned_seqs'.")
    if 'prediction' not in df.columns:
        df[f'prediction_{task}'] = None

    for i, row in df.iterrows():
        if isinstance(row[f'prediction_{task}'], str):
            continue
        try:
            seqs = literal_eval(row['unaligned_seqs'])
            prompt = prepare_input_string(seqs, SPECIAL_TOKENS[task])
            result = generate_output(model, tokenizer, prompt, task)
            print(f"Result: {result}")
            df.at[i, f'prediction_{task}'] = result
            print(f"Processed row {i+1}/{len(df)}")
        except Exception as e:
            print(f"Error processing row {i}: {e}")

    df.to_csv(output_path, index=False)
    print(f"Saved predictions to: {output_path}")

def handle_fasta_input(model, tokenizer, path, output_path, task):
    print(f"Processing FASTA: {path}")
    seqs = [seq for _, seq in fasta_iter(path)]
    prompt = prepare_input_string(seqs, SPECIAL_TOKENS[task])
    result = generate_output(model, tokenizer, prompt, task)
    print(f"Result: {result}")
    with open(output_path, 'w') as f:
        f.write(">prediction\n")
        f.write(result + "\n")
    print(f"Saved prediction to: {output_path}")

def handle_fasta_folder(model, tokenizer, folder_path, output_folder, task):
    print(f"Processing folder of FASTA files: {folder_path}")
    os.makedirs(output_folder, exist_ok=True)
    fasta_files = glob.glob(os.path.join(folder_path, "*.fasta"))
    if not fasta_files:
        raise FileNotFoundError("No .fasta files found in the specified folder.")
    for fasta_path in fasta_files:
        try:
            seqs = [seq for _, seq in fasta_iter(fasta_path)]
            prompt = prepare_input_string(seqs, SPECIAL_TOKENS[task])
            result = generate_output(model, tokenizer, prompt, task)
            print(f"Result: {result}")
            out_path = os.path.join(output_folder, os.path.basename(fasta_path))
            with open(out_path, 'w') as f:
                f.write(">prediction\n")
                f.write(result + "\n")
            print(f"Saved: {out_path}")
        except Exception as e:
            print(f"Failed on {fasta_path}: {e}")

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

    model, tokenizer = load_model_and_tokenizer(args.model_path)

    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input path does not exist: {args.input}")

    if os.path.isdir(args.input):
        handle_fasta_folder(model, tokenizer, args.input, args.output, task)
    elif args.input.endswith(".csv"):
        handle_csv_input(model, tokenizer, args.input, args.output, task)
    elif args.input.endswith((".fasta", ".fa")):
        handle_fasta_input(model, tokenizer, args.input, args.output, task)
    else:
        raise ValueError("Unsupported input type. Use a .csv, .fasta file, or folder of .fasta files.")

if __name__ == "__main__":
    main(parse_arguments())