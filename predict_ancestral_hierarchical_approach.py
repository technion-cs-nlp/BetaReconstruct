import os
import json
import random
import math
import argparse
from itertools import groupby

import torch
from typing import List, Dict, Tuple, Optional
from transformers import ZambaForCausalLM, GPT2TokenizerFast, AutoConfig

# Outgroup clades (marsupials), as in the fine-tuning data
DEFAULT_OUTGROUP_CLADES = ["Dasyuromorphia", "Didelphimorphia", "Diprotodontia", "Peramelemorphia", "Dromiciops"]
OUTGROUP_COUNT = 2   # the model was fine-tuned with exactly 2 outgroup sequences
MAX_INGROUP = 8      # the model was fine-tuned on exactly 8 ingroup sequences
MIN_CLADE_SIZE = 5   # clades with fewer sequences are merged into the "other" clade
OTHER_CLADE = "other"
ROOT_CLADE = "root"

def fasta_iter(fasta_name: str) -> Dict[str, str]:
    """Parse FASTA file into a dict of header -> sequence."""
    with open(fasta_name) as fh:
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        return {next(header)[1:].strip(): "".join(s.strip() for s in next(faiter)) for header in faiter}

def normalize_name(name: str) -> str:
    return "_".join(name.split()).lower()

def load_species_list(species_list_folder: str, outgroup_clades: List[str]) -> Tuple[Dict[str, List[str]], set]:
    """
    Load species list files (<clade>.txt, one species per line) and normalize names.
    Species of the outgroup clades are kept apart: they are never reconstructed, only used as outgroup.
    """
    species_dict = {OTHER_CLADE: []}
    outgroup_species = set()
    found_outgroup_clades = set()
    for root, _, files in os.walk(species_list_folder):
        for fname in sorted(files):
            if not fname.endswith(".txt"):
                continue
            with open(os.path.join(root, fname)) as f:
                species = [normalize_name(line) for line in f if line.strip()]
            clade = fname.split(".")[0]
            if clade in outgroup_clades:
                outgroup_species.update(species)
                found_outgroup_clades.add(clade)
            else:
                species_dict[clade] = species
    missing = [clade for clade in outgroup_clades if clade not in found_outgroup_clades]
    if missing:
        print(f"Warning: no species list for the outgroup clades {missing}; their species will be treated as ingroup.")
    return species_dict, outgroup_species

def assign_to_clades(records: List[Tuple[str, str]], species_map: Dict[str, List[str]]) -> Dict[str, List[Tuple[str, str]]]:
    """Assign (name, sequence) records to the clades of species_map; unmatched species go to the "other" clade."""
    clades_data = {clade: [] for clade in species_map}
    for key, seq in records:
        norm_key = key.lower()
        for clade, species_list in species_map.items():
            if norm_key in species_list:
                clades_data[clade].append((key, seq))
                break
        else:
            clades_data[OTHER_CLADE].append((key, seq))
    return clades_data

def merge_small_clades(clades_data: Dict[str, List[Tuple[str, str]]], min_clade_size: int) -> Dict[str, List[Tuple[str, str]]]:
    """Merge clades with fewer than min_clade_size sequences into the "other" clade."""
    merged = {OTHER_CLADE: list(clades_data.get(OTHER_CLADE, []))}
    for clade, data in clades_data.items():
        if clade == OTHER_CLADE:
            continue
        if len(data) < min_clade_size:
            merged[OTHER_CLADE].extend(data)
        else:
            merged[clade] = list(data)
    return {clade: data for clade, data in merged.items() if data}

def divide_species_to_groups(clades_data: Dict[str, List[Tuple[str, str]]], level: int, max_group_size: int, min_clade_size: int) -> Dict[str, List[Tuple[str, str]]]:
    """
    Balanced split capped at max_group_size: a clade with n sequences is split into
    ceil(n / max_group_size) random groups whose sizes differ by at most one
    (e.g., 20 -> 7+7+6, 17 -> 6+6+5). Clades with fewer than min_clade_size sequences form a single group.
    """
    reconstruction_groups = {}
    for clade, data in clades_data.items():
        if not data:
            continue
        data = list(data)
        random.shuffle(data)
        num_groups = 1 if len(data) < min_clade_size else math.ceil(len(data) / max_group_size)
        base_size, remainder = divmod(len(data), num_groups)
        start = 0
        for group_index in range(num_groups):
            size = base_size + (1 if group_index < remainder else 0)
            reconstruction_groups[f"LVL_{level}_clade_{clade}_group_{group_index}"] = data[start:start + size]
            start += size
    return reconstruction_groups

def reconstruct_internal_node(model, tokenizer, device, ingroup_seqs: List[str], outgroup_pool: List[Tuple[str, str]], outgroup_count: int) -> Tuple[str, List[str]]:
    """
    Reconstruct the ancestor of the ingroup sequences, given outgroup_count outgroup sequences
    sampled from the outgroup pool. Input format is identical to training:
        out_1#out_2#in_1|in_2|...|in_n<RECONSTRUCT>
    Returns the prediction and the names of the outgroup sequences that were used.
    """
    selected_outgroups = random.sample(outgroup_pool, outgroup_count)
    outgroup_seqs = [seq for _, seq in selected_outgroups]
    input_str = "".join(seq + "#" for seq in outgroup_seqs) + "|".join(ingroup_seqs) + "<RECONSTRUCT>"

    encoded = tokenizer(input_str, return_tensors='pt')
    encoded = {k: v.to(device) for k, v in encoded.items()}

    # Stop generation at EOS or at the first [UNK] token
    stop_token_ids = [tokenizer.eos_token_id]
    if tokenizer.unk_token_id is not None and tokenizer.unk_token_id != tokenizer.eos_token_id:
        stop_token_ids.append(tokenizer.unk_token_id)

    with torch.no_grad():
        try:
            output = model.generate(**encoded, tokenizer=tokenizer, return_dict_in_generate=True, max_length=model.config.max_position_embeddings, eos_token_id=stop_token_ids)
        except Exception:
            output = model.generate(**encoded, tokenizer=tokenizer, return_dict_in_generate=True, max_new_tokens=650, eos_token_id=stop_token_ids)

    decoded = tokenizer.decode(output.sequences[0])
    prediction = decoded.split("<RECONSTRUCT>")[1]
    # Keep the protein sequence only up to the first [UNK]
    if tokenizer.unk_token is not None and tokenizer.unk_token in prediction:
        prediction = prediction.split(tokenizer.unk_token)[0]
    prediction = "".join(prediction.split())
    prediction = prediction.split("|")[0].split("#")[0]
    prediction = prediction.replace(tokenizer.decode(tokenizer.eos_token_id), "")
    return prediction, [name for name, _ in selected_outgroups]

def reconstruct_hierarchically(model, tokenizer, ingroup_by_clade: Dict[str, List[Tuple[str, str]]], outgroup_pool: List[Tuple[str, str]], args) -> Dict:
    """
    Level by level: the sequences of each clade are split into groups of up to max_group_size,
    the ancestor of each group is reconstructed, and these ancestors are the input of the next level.
    When all remaining sequences fit in one group, its ancestor (the ancestor of the whole ingroup) is the prediction.
    """
    results = {}
    outgroups_used = {}
    current_data = ingroup_by_clade
    level = 0
    while True:
        current_data = merge_small_clades(current_data, args.min_clade_size)
        num_seqs = sum(len(data) for data in current_data.values())
        is_last_level = num_seqs <= args.max_group_size
        if is_last_level:
            groups = {f"LVL_{level}_clade_{ROOT_CLADE}_group_0": [item for data in current_data.values() for item in data]}
        else:
            groups = divide_species_to_groups(current_data, level, args.max_group_size, args.min_clade_size)
        results[f"reconstruction_groups_{level}"] = groups
        print(f"Level {level}: {num_seqs} sequences in {len(groups)} groups")

        next_data = {}
        for title, group in groups.items():
            clade = title.split(f"LVL_{level}_clade_")[1].rsplit("_group_", 1)[0]
            keys = [k for k, _ in group]
            seqs = [s for _, s in group]
            prediction, outgroups_used[title] = reconstruct_internal_node(model, tokenizer, model.device, seqs, outgroup_pool, args.outgroup_count)
            next_data.setdefault(clade, []).append(("-".join(keys), prediction))

        if is_last_level:
            return {**results, "outgroup_pool": [k for k, _ in outgroup_pool], "outgroups_used": outgroups_used, "prediction": prediction}
        current_data = next_data
        level += 1

def find_families(fasta_files_folder: str) -> List[Tuple[str, str, Optional[str]]]:
    """
    Families to reconstruct, as (family_name, ingroup_or_all_fasta, outgroup_fasta):
    - a subfolder with ingroup.fasta and outgroup.fasta (as in data/mammalian_validation): the outgroup is given;
    - any other .fasta file: the outgroup is identified by the species lists of the outgroup clades.
    """
    families = []
    for path, _, files in os.walk(fasta_files_folder):
        if "ingroup.fasta" in files and "outgroup.fasta" in files:
            families.append((os.path.basename(os.path.normpath(path)), os.path.join(path, "ingroup.fasta"), os.path.join(path, "outgroup.fasta")))
            files = [f for f in files if f not in ("ingroup.fasta", "outgroup.fasta")]
        for fasta_file in sorted(files):
            if fasta_file.endswith(".fasta"):
                families.append((fasta_file.split(".")[0], os.path.join(path, fasta_file), None))
    return sorted(families)

def read_family(fasta_path: str, outgroup_path: Optional[str], outgroup_species: set) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]]]:
    """Return the ingroup and outgroup (name, unaligned sequence) records of a family."""
    records = [(key, seq.replace("-", "")) for key, seq in fasta_iter(fasta_path).items()]
    if outgroup_path is not None:
        outgroup = [(key, seq.replace("-", "")) for key, seq in fasta_iter(outgroup_path).items()]
        return records, outgroup
    ingroup = [(key, seq) for key, seq in records if key.lower() not in outgroup_species]
    outgroup = [(key, seq) for key, seq in records if key.lower() in outgroup_species]
    return ingroup, outgroup

def initialize_model(args):
    """Load model and tokenizer based on CLI args."""
    tokenizer = GPT2TokenizerFast.from_pretrained(args.model_path)

    try:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        if device.type == "cuda":
            print("Using GPU")
            model = ZambaForCausalLM.from_pretrained(args.model_path)
        else:
            print("Using CPU")
            config = AutoConfig.from_pretrained(args.model_path)
            config.use_mamba_kernels = False
            model = ZambaForCausalLM.from_pretrained(args.model_path, config=config)
    except Exception as e:
        exit(f"Failed to load model - {e}, exit")

    model.generation_config.pad_token_id = tokenizer.pad_token_id
    model.generation_config.eos_token_id = tokenizer.eos_token_id
    model = model.to(device)
    model.eval()
    return model, tokenizer

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model_path", help="model path or name to load", type=str, required=True)
    parser.add_argument("--fasta_files_folder", help="folder with the families to reconstruct: .fasta files (one per family), "
                        "or subfolders with ingroup.fasta and outgroup.fasta (one per family). Species names must match the names in 'species_list_folder'", type=str, required=True)
    parser.add_argument("--results_folder", help="path to save results", type=str, required=True)
    parser.add_argument("--species_list_folder", help="folder containing txt files, where each txt file contains the species of a specific clade (e.g., an order)", type=str, required=True)
    parser.add_argument("--outgroup_clades", help="comma-separated names of the clades (species list files) that form the outgroup; "
                        "used for .fasta files without a separate outgroup.fasta", default=",".join(DEFAULT_OUTGROUP_CLADES), type=str)
    parser.add_argument("--outgroup_count", help="number of outgroup sequences in each model input", default=OUTGROUP_COUNT, type=int)
    parser.add_argument("--max_group_size", help="maximal number of ingroup sequences in each model input", default=MAX_INGROUP, type=int)
    parser.add_argument("--min_clade_size", help="clades with fewer sequences are merged into the 'other' clade", default=MIN_CLADE_SIZE, type=int)
    parser.add_argument("--seed", help="random seed (grouping and outgroup sampling)", default=42, type=int)
    parser.add_argument("--model_name", help="name to save results with", default="", type=str)
    return parser.parse_args()

def main(args):
    print(args)
    random.seed(args.seed)
    os.makedirs(args.results_folder, exist_ok=True)
    outgroup_clades = [clade.strip() for clade in args.outgroup_clades.split(",") if clade.strip()]
    species_map, outgroup_species = load_species_list(args.species_list_folder, outgroup_clades)
    families = find_families(args.fasta_files_folder)
    if not families:
        raise FileNotFoundError("No .fasta files (or subfolders with ingroup.fasta and outgroup.fasta) found in the specified folder.")
    model, tokenizer = initialize_model(args)

    for family_name, fasta_path, outgroup_path in families:
        print(f"Processing: {family_name}")
        ingroup, outgroup_pool = read_family(fasta_path, outgroup_path, outgroup_species)
        print(f"{len(ingroup)} ingroup and {len(outgroup_pool)} outgroup sequences: {[k for k, _ in outgroup_pool]}")
        if len(outgroup_pool) < args.outgroup_count:
            print(f"[SKIP] {family_name}: fewer than {args.outgroup_count} outgroup sequences.")
            continue
        if not ingroup:
            print(f"[SKIP] {family_name}: no ingroup sequences.")
            continue

        ingroup_by_clade = assign_to_clades(ingroup, species_map)
        results = reconstruct_hierarchically(model, tokenizer, ingroup_by_clade, outgroup_pool, args)

        if args.model_name != "":
            output_path = os.path.join(args.results_folder, f'{family_name}_{args.model_name}.json')
        else:
            output_path = os.path.join(args.results_folder, f'{family_name}.json')
        with open(output_path, 'w') as f:
            json.dump(results, f)
        print(f"Saved: {output_path}")

if __name__ == "__main__":
    main(parse_arguments())
