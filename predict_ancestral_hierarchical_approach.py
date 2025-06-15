import os
import json
import random
import math
import argparse
from itertools import groupby

import torch
import pandas as pd
from typing import List, Dict, Tuple
from transformers import ZambaForCausalLM, GPT2TokenizerFast

def fasta_iter(fasta_name: str) -> Dict[str, str]:
    """Parse FASTA file into a dict of header -> sequence."""
    with open(fasta_name) as fh:
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        return {next(header)[1:].strip(): "".join(s.strip() for s in next(faiter)) for header in faiter}

def divide_species_to_groups(orders_and_data: Dict[str, List[Tuple[str, str]]], level: int) -> Dict[str, List[Tuple[str, str]]]:
    """Group species data into reconstruction groups by size."""
    reconstruction_groups = {}
    # count number of seqs in input
    num_seqs_in_dictionary = 0
    for clade, data in orders_and_data.items():
        num_seqs_in_dictionary += len(data)
    print(f'num_seqs_in_dictionary = {num_seqs_in_dictionary}')
    other_dictionary = {}
    # merge all seqs in order_and_data to other if entire number below 10
    if num_seqs_in_dictionary < 10:
        for clade, data in orders_and_data.items():
            seqs = other_dictionary.get(f"LVL_{level}_clade_other_group_0", [])
            seqs.extend(data)
            other_dictionary[f"LVL_{level}_clade_other_group_0"] = seqs
            other_dictionary[clade] = []
        orders_and_data = other_dictionary
    
    # merge clades to groups
    for clade, data in orders_and_data.items():
        group_data = data.copy()
        for num_species in [9, 8, 7, 6, 10, 5]:
            if len(group_data) % num_species == 0 or len(group_data) % num_species >= 5:
                num_groups = 0
                while group_data:
                    group = [group_data.pop(random.randrange(len(group_data))) for _ in range(min(num_species, len(group_data)))]
                    group_key = f"LVL_{level}_clade_{clade}_group_{num_groups}"
                    reconstruction_groups[group_key] = group
                    num_groups += 1
                break
            elif len(group_data) < 5:
                reconstruction_groups[f"LVL_{level}_clade_{clade}_group_0"] = group_data
                break
    return reconstruction_groups

def run_model_on_group(model, tokenizer, device, group: List[Tuple[str, str]], clade: str) -> Tuple[str, str]:
    """Generate reconstruction prediction for a group."""
    keys = [k for k, _ in group]
    seqs = [s for _, s in group]
    input_str = "|".join(seqs) + "<RECONSTRUCT>"

    encoded = tokenizer(input_str, return_tensors='pt')
    encoded = {k: v.to(device) for k, v in encoded.items()}

    try:
        output = model.generate(**encoded, tokenizer=tokenizer, return_dict_in_generate=True, max_length=model.config.max_position_embeddings)
    except Exception:
        output = model.generate(**encoded, tokenizer=tokenizer, return_dict_in_generate=True, max_new_tokens=650)
    
    decoded = tokenizer.decode(output.sequences[0]).replace(" ", "")
    prediction = decoded.split("<RECONSTRUCT>")[1].split("|")[0].replace(tokenizer.decode(tokenizer.eos_token_id), "").replace(tokenizer.decode(tokenizer.unk_token_id), "")
    print(f'x = {input_str}')
    print()
    print(f'p = {prediction}')
    print()
    return "-".join(keys), prediction

def load_species_list(species_list_folder: str) -> Dict[str, List[str]]:
    """Load species list files and normalize names."""
    species_dict = {"other": []}
    for root, _, files in os.walk(species_list_folder):
        for fname in files:
            if fname.endswith(".txt"):
                with open(os.path.join(root, fname)) as f:
                    species = ["_".join(line.split()).lower() for line in f]
                species_dict[fname.split(".")[0]] = species
    return species_dict

def initialize_model(args):
    """Load model and tokenizer based on CLI args."""
    tokenizer = GPT2TokenizerFast.from_pretrained(args.model_path)

    try:
        model = ZambaForCausalLM.from_pretrained(args.model_path)
    except Exception as e:
        exit(f"Failed to load model - {e}, exit")

    return model.to("cuda" if torch.cuda.is_available() else "cpu"), tokenizer

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model_path", help="model path or name to load", type=str)
    parser.add_argument("--fasta_files_folder", help="folder with fasta files to predict ancestral. make sure the specie names matches the name of the 'species_list_folder'", type=str)
    parser.add_argument("--results_folder", help="path to save results", type=str)
    parser.add_argument("--species_list_folder", help="folder containing txt files, where each txt file contains the species of a specific order", type=str)
    parser.add_argument("--model_name", help="name to save results with", default="", type=str)
    return parser.parse_args()

def main(args):
    print(args)
    os.makedirs(args.results_folder, exist_ok=True)
    model, tokenizer = initialize_model(args)
    print(tokenizer)
    print(model)
    species_map = load_species_list(args.species_list_folder)

    for path, _, files in os.walk(args.fasta_files_folder):
        for fasta_file in files:
            if not fasta_file.endswith(".fasta"):
                continue

            family_name = fasta_file.split(".")[0]
            print(f"Processing: {fasta_file}")
            data = fasta_iter(os.path.join(path, fasta_file))
            orders_data = {order: [] for order in species_map}

            for key, seq in list(data.items()):
                norm_key = key.lower()
                found = False
                for family, species_list in species_map.items():
                    if norm_key in species_list:
                        orders_data[family].append((key, seq.replace("-", "")))
                        found = True
                        break
                if not found:
                    orders_data["other"].append((key, seq.replace("-", "")))

            for k in list(orders_data):
                if len(orders_data[k]) < 5 and k != "other":
                    orders_data["other"].extend(orders_data[k])
                    orders_data[k] = []

            # Dynamic level reconstruction
            current_level = 0
            all_reconstruction_groups = {}
            current_data = orders_data
            prediction = None
            
            while True:
                print(f'current_level = {current_level}')
                groups = divide_species_to_groups(current_data, current_level)
                print(f'len(groups) = {len(groups)}')
                #print(f'current_data = {current_data}')
                all_reconstruction_groups[f"reconstruction_groups_{current_level}"] = groups

                if len(groups) == 1:
                    # Root ancestor is found
                    only_group = next(iter(groups.values()))
                    prediction = run_model_on_group(model, tokenizer, model.device, only_group, clade="root")[1]
                    break

                next_data = {}
                for title, group in groups.items():
                    clade = title.split('_')[3]
                    group_id, pred = run_model_on_group(model, tokenizer, model.device, group, clade)
                    next_data.setdefault(clade, []).append((group_id, pred))
                
                current_data = next_data
                current_level += 1
            if args.model_name != "":
                output_path = os.path.join(args.results_folder, f'{family_name}_{args.model_name}.json')
            else:
                output_path = os.path.join(args.results_folder, f'{family_name}.json')
            with open(output_path, 'w') as f:
                json.dump({
                    **all_reconstruction_groups,
                    "prediction": prediction
                }, f)

if __name__ == "__main__":
    main(parse_arguments())