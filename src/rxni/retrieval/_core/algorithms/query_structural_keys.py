"""Frozen Library-v4 query structural representations.

Behavior-preserving extraction from:

    template_library_v4_builder_v0_7_1/
        build_template_library_v4.py
        build_membership_v4.py

The chemistry and hashing behavior in this module must remain identical to the
builder that generated the v4 membership library.

Do not optimize or chemically modify these functions without releasing a new
representation version and validating full query/member parity.
"""

from __future__ import annotations

import hashlib
import math
import re

from collections import Counter, defaultdict
from typing import Any, Iterable

def h64(text: str) -> int:
    return int.from_bytes(hashlib.blake2b(text.encode("utf-8"), digest_size=8).digest(), "little", signed=False)


def require_rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem import rdChemReactions
    except Exception as exc:
        raise RuntimeError(
            "RDKit is required for requested chemistry-derived features. "
            "Run in the aizynthfinder/RDKit environment."
        ) from exc
    return Chem, rdChemReactions


def bond_order_label(bond) -> str:
    try:
        v = float(bond.GetBondTypeAsDouble())
        return f"{v:g}"
    except Exception:
        return str(bond.GetBondType())


def atom_label(atom) -> str:
    z = atom.GetAtomicNum()
    sym = atom.GetSymbol() if z else "*"
    aromatic = "a" if atom.GetIsAromatic() else ""
    charge = atom.GetFormalCharge()
    ch = f"{charge:+d}" if charge else ""
    return f"{sym}{aromatic}{ch}"


def side_maps(mols):
    out = {}
    comp_by_map = {}
    atom_by_map = {}
    for ci, mol in enumerate(mols):
        for a in mol.GetAtoms():
            m = a.GetAtomMapNum()
            if m:
                out[m] = atom_label(a)
                comp_by_map[m] = ci
                atom_by_map[m] = a
    return out, comp_by_map, atom_by_map


def side_bonds(mols):
    d = {}
    bond_objects = {}
    for ci, mol in enumerate(mols):
        for b in mol.GetBonds():
            a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
            m1, m2 = a1.GetAtomMapNum(), a2.GetAtomMapNum()
            if not m1 or not m2:
                continue
            key = tuple(sorted((m1, m2)))
            d[key] = bond_order_label(b)
            bond_objects[key] = (ci, b)
    return d, bond_objects


def graph_cycle_rank(mols) -> int:
    # Sum(E - V + C) over disconnected side components.
    # Template components are provided separately.
    total = 0
    for mol in mols:
        v = mol.GetNumAtoms()
        e = mol.GetNumBonds()
        # GetMolFrags handles any internal disconnected fragments safely.
        try:
            from rdkit import Chem
            c = len(Chem.GetMolFrags(mol))
        except Exception:
            c = 1
        total += max(0, e - v + c)
    return total


def _wl_signature(mol, node_seed: dict[int, str], keep: set[int] | None = None, rounds: int = 3) -> str:
    """Compact map-number-invariant local graph signature.

    Atom mapping is used upstream to locate the reaction center, but map IDs are
    intentionally excluded from the signature itself.
    """
    if keep is None:
        keep = set(range(mol.GetNumAtoms()))
    labels = {i: h64(node_seed[i]) for i in keep}
    for _ in range(rounds):
        new = {}
        for i in keep:
            a = mol.GetAtomWithIdx(i)
            nbr = []
            for b in a.GetBonds():
                j = b.GetOtherAtomIdx(i)
                if j in keep:
                    nbr.append((bond_order_label(b), labels[j]))
            nbr.sort()
            payload = f"{labels[i]}|" + "|".join(f"{bo}:{lab}" for bo, lab in nbr)
            new[i] = h64(payload)
        labels = new

    node_part = ",".join(str(x) for x in sorted(labels.values()))
    edge_part = []
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in keep and j in keep:
            x, y = sorted((labels[i], labels[j]))
            edge_part.append(f"{x}-{bond_order_label(b)}-{y}")
    edge_part.sort()
    return f"N[{node_part}]E[{','.join(edge_part)}]"


def component_signature(mol) -> str:
    seed = {
        a.GetIdx(): f"{atom_label(a)}|deg={a.GetDegree()}"
        for a in mol.GetAtoms()
    }
    return _wl_signature(mol, seed, rounds=3)


def shortest_path_len(mol, idx1: int, idx2: int):
    try:
        from rdkit import Chem
        path = Chem.GetShortestPath(mol, idx1, idx2)
        return len(path) - 1 if path else None
    except Exception:
        return None


def atom_idx_for_map(mol, amap):
    for a in mol.GetAtoms():
        if a.GetAtomMapNum() == amap:
            return a.GetIdx()
    return None


def local_center_signature(mols, changed_maps: set[int], radius: int) -> str:
    pieces = []
    for ci, mol in enumerate(mols):
        dist = {}
        frontier = []
        for a in mol.GetAtoms():
            if a.GetAtomMapNum() in changed_maps:
                dist[a.GetIdx()] = 0
                frontier.append(a.GetIdx())

        for d in range(radius):
            nxt = []
            for idx in frontier:
                a = mol.GetAtomWithIdx(idx)
                for nb in a.GetNeighbors():
                    j = nb.GetIdx()
                    if j not in dist:
                        dist[j] = d + 1
                        nxt.append(j)
            frontier = nxt

        keep = set(dist)
        if not keep:
            pieces.append(f"C{ci}:EMPTY")
            continue

        # Mapping IDs are deliberately absent. Distance from the mapped edit
        # center plus atom/bond topology defines the local structural signature.
        seed = {}
        for idx in keep:
            a = mol.GetAtomWithIdx(idx)
            seed[idx] = (
                f"{atom_label(a)}|d={dist[idx]}|"
                f"center={int(dist[idx] == 0)}|deg={a.GetDegree()}"
            )

        pieces.append(f"C{ci}:{_wl_signature(mol, seed, keep=keep, rounds=3)}")

    # Components are unordered for signature purposes.
    return "||".join(sorted(pieces))


def derive_features_from_mols(reactants, products, store_readable=False) -> dict[str, Any]:
    rmap, rcomp, ratoms = side_maps(reactants)
    pmap, pcomp, patoms = side_maps(products)
    rbonds, _ = side_bonds(reactants)
    pbonds, _ = side_bonds(products)

    deleted = []
    added = []
    changed = []

    for k, bo in rbonds.items():
        if k not in pbonds:
            deleted.append((k, bo))
        elif pbonds[k] != bo:
            changed.append((k, bo, pbonds[k]))
    for k, bo in pbonds.items():
        if k not in rbonds:
            added.append((k, bo))

    changed_maps = set()
    for (a, b), _ in deleted + added:
        changed_maps.update((a, b))
    for (a, b), _, _ in changed:
        changed_maps.update((a, b))

    def pair_labels(pair, amap):
        a, b = pair
        return tuple(sorted((amap.get(a, "?"), amap.get(b, "?"))))

    ed_tokens = []
    broken_elements = []
    formed_elements = []
    changed_elements = []

    for pair, bo in deleted:
        labs = pair_labels(pair, rmap)
        ed_tokens.append(f"DEL:{labs[0]}-{labs[1]}:{bo}")
        broken_elements.append(f"{labs[0]}-{labs[1]}")
    for pair, bo in added:
        labs = pair_labels(pair, pmap)
        ed_tokens.append(f"ADD:{labs[0]}-{labs[1]}:{bo}")
        formed_elements.append(f"{labs[0]}-{labs[1]}")
    for pair, old, new in changed:
        labs = pair_labels(pair, rmap | pmap)
        ed_tokens.append(f"CHG:{labs[0]}-{labs[1]}:{old}>{new}")
        changed_elements.append(f"{labs[0]}-{labs[1]}:{old}>{new}")
    ed_tokens.sort()
    edit_sig = "|".join(ed_tokens) if ed_tokens else "NO_MAPPED_BOND_EDIT"

    # Determine participating reactant/product components.
    r_part = {rcomp[m] for m in changed_maps if m in rcomp}
    p_part = {pcomp[m] for m in changed_maps if m in pcomp}

    intermolecular = False
    intramolecular = False
    new_ring_sizes = []
    for pair, _bo in added:
        a, b = pair
        if a in rcomp and b in rcomp:
            if rcomp[a] != rcomp[b]:
                intermolecular = True
            else:
                intramolecular = True
                mol = reactants[rcomp[a]]
                ia = atom_idx_for_map(mol, a)
                ib = atom_idx_for_map(mol, b)
                if ia is not None and ib is not None:
                    d = shortest_path_len(mol, ia, ib)
                    if d is not None and d >= 1:
                        new_ring_sizes.append(d + 1)

    ring_formed = bool(new_ring_sizes)
    # Ring-broken proxy: deleted mapped bond whose endpoints remain in one product component.
    ring_broken = False
    for pair, _bo in deleted:
        a, b = pair
        if a in pcomp and b in pcomp and pcomp[a] == pcomp[b]:
            # If endpoints are still connected through another path on product side,
            # deleting the original bond opened/rewired a cycle.
            mol = products[pcomp[a]]
            ia = atom_idx_for_map(mol, a)
            ib = atom_idx_for_map(mol, b)
            if ia is not None and ib is not None:
                d = shortest_path_len(mol, ia, ib)
                if d is not None and d >= 1:
                    ring_broken = True

    rcomp_sigs = [component_signature(m) for m in reactants]
    pcomp_sigs = [component_signature(m) for m in products]
    rcomp_hashes = [h64(x) for x in rcomp_sigs]
    pcomp_hashes = [h64(x) for x in pcomp_sigs]

    topo_sig = (
        f"R{len(reactants)}P{len(products)}|"
        f"RPART{len(r_part)}|PPART{len(p_part)}|"
        f"INTRA{int(intramolecular)}|INTER{int(intermolecular)}|"
        f"MERGE{max(0, len(r_part)-len(p_part))}|"
        f"SPLIT{max(0, len(p_part)-len(r_part))}"
    )
    # Direction-sensitive, multiplicity-preserving component evidence.  Unlike
    # component_topology_hash, this retains the identities of all components and
    # is therefore suitable as the compact System-I component representation.
    component_evidence_sig = (
        "R[" + ",".join(str(x) for x in sorted(rcomp_hashes)) + "]>>P["
        + ",".join(str(x) for x in sorted(pcomp_hashes)) + "]"
    )

    raw = {}
    for radius in (0, 1, 2):
        sig = local_center_signature(reactants, changed_maps, radius) + ">>" + local_center_signature(products, changed_maps, radius)
        raw[radius] = sig

    # Conservative edit-driven FGI label: useful as a bucket, not semantic truth.
    fgi = infer_fgi_from_edit_tokens(ed_tokens)

    out = {
        "n_reactant_atoms": sum(m.GetNumAtoms() for m in reactants),
        "n_product_atoms": sum(m.GetNumAtoms() for m in products),
        "n_mapped_atoms": len(set(rmap) & set(pmap)),
        "n_changed_atoms": len(changed_maps),
        "n_reactant_components": len(reactants),
        "n_product_components": len(products),
        "n_bond_additions": len(added),
        "n_bond_deletions": len(deleted),
        "n_bond_order_changes": len(changed),
        "atom_conservation_fraction": (
            len(set(rmap) & set(pmap)) / max(len(set(rmap) | set(pmap)), 1)
        ),
        "component_delta": len(products) - len(reactants),
        "participating_component_count": len(r_part),
        "intramolecular_edit": bool(intramolecular),
        "intermolecular_edit": bool(intermolecular),
        "component_merge_count": max(0, len(r_part) - len(p_part)),
        "component_split_count": max(0, len(p_part) - len(r_part)),
        "edit_signature_hash": h64(edit_sig),
        "bond_transformation_signature": edit_sig,
        "formed_bond_elements": "|".join(sorted(formed_elements)),
        "broken_bond_elements": "|".join(sorted(broken_elements)),
        "changed_bond_elements": "|".join(sorted(changed_elements)),
        "raw_center_hash_r0": h64(raw[0]),
        "raw_center_hash_r1": h64(raw[1]),
        "raw_center_hash_r2": h64(raw[2]),
        "reactant_component_hashes": rcomp_hashes,
        "product_component_hashes": pcomp_hashes,
        "component_topology_hash": h64(topo_sig),
        "component_evidence_hash": h64(component_evidence_sig),
        "ring_formed": ring_formed,
        "ring_broken": ring_broken,
        "ring_count_delta": graph_cycle_rank(products) - graph_cycle_rank(reactants),
        "new_ring_sizes": sorted(set(new_ring_sizes)),
        "fgi_signature": fgi,
        "fgi_signature_hash": h64(fgi),
        "_edit_tokens": ed_tokens,
        "_raw_tokens": [raw[0], raw[1], raw[2]],
        "_component_topology_signature": topo_sig,
    }

    if store_readable:
        out["edit_signature"] = edit_sig
        out["raw_center_signature_r0"] = raw[0]
        out["raw_center_signature_r1"] = raw[1]
        out["raw_center_signature_r2"] = raw[2]
        out["reactant_component_signatures"] = rcomp_sigs
        out["product_component_signatures"] = pcomp_sigs
        out["component_topology_signature"] = topo_sig
        out["component_evidence_signature"] = component_evidence_sig

    return out


def infer_fgi_from_edit_tokens(tokens: list[str]) -> str:
    # Intentionally conservative pilot buckets.
    t = "|".join(tokens)
    rules = [
        (("DEL:Br-C", "ADD:C-N"), "c_br_to_n"),
        (("DEL:C-Br", "ADD:C-N"), "c_br_to_n"),
        (("DEL:C-Cl", "ADD:C-N"), "c_cl_to_n"),
        (("DEL:C-F", "ADD:C-N"), "c_f_to_n"),
        (("DEL:C-I", "ADD:C-N"), "c_i_to_n"),
        (("DEL:C-O", "ADD:C-N"), "c_o_to_n"),
        (("DEL:C-Br", "ADD:C-O"), "c_br_to_o"),
        (("DEL:C-Cl", "ADD:C-O"), "c_cl_to_o"),
        (("DEL:C-O", "ADD:C-C"), "c_o_to_c"),
        (("DEL:Br-C", "ADD:C-C"), "c_br_to_c"),
        (("DEL:C-Br", "ADD:C-C"), "c_br_to_c"),
        (("DEL:C-Cl", "ADD:C-C"), "c_cl_to_c"),
    ]
    for required, label in rules:
        if all(x in t for x in required):
            return label
    # Bond-order signatures retain useful directional information.
    chg = [x for x in tokens if x.startswith("CHG:")]
    if chg:
        return "bond_order_change:" + ";".join(chg)
    if tokens:
        return "edit:" + ";".join(tokens)
    return "none"


def _parse_side(Chem, text: str):
    mols = []
    for comp in (text or "").split("."):
        comp = comp.strip()
        if not comp:
            continue
        mol = Chem.MolFromSmiles(comp)
        if mol is None:
            mol = Chem.MolFromSmiles(comp, sanitize=False)
        if mol is None:
            raise ValueError(f"Could not parse reaction component: {comp[:120]!r}")
        mols.append(mol)
    return mols


def parse_mapped_reaction(reaction: str):
    Chem, _ = require_rdkit()
    text = str(reaction).strip()
    if text.count(">") == 2:
        lhs, _agents, rhs = text.split(">", 2)
    elif ">>" in text:
        lhs, rhs = text.split(">>", 1)
    else:
        raise ValueError("Reaction must contain '>>' or reactants>agents>products")
    reactants = _parse_side(Chem, lhs)
    products = _parse_side(Chem, rhs)
    if not reactants or not products:
        raise ValueError("Reaction has empty reactant or product side")
    return reactants, products



def derive_query_structural_keys(
    mapped_reaction: str,
    *,
    store_readable: bool = False,
) -> dict[str, Any]:
    """Derive the concrete forward structural representation used by Library v4."""

    reactants, products = parse_mapped_reaction(
        mapped_reaction
    )

    f = derive_features_from_mols(
        reactants,
        products,
        store_readable=store_readable,
    )

    out = {
        "forward_raw_center_hash_r0":
            f["raw_center_hash_r0"],

        "forward_raw_center_hash_r1":
            f["raw_center_hash_r1"],

        "forward_raw_center_hash_r2":
            f["raw_center_hash_r2"],

        "forward_edit_signature_hash":
            f["edit_signature_hash"],

        "forward_component_evidence_hash":
            f["component_evidence_hash"],

        "forward_component_topology_hash":
            f.get("component_topology_hash"),

        "forward_reactant_component_hashes":
            f.get("reactant_component_hashes"),

        "forward_product_component_hashes":
            f.get("product_component_hashes"),
    }

    # Canonical System-I v4 member aliases.
    out["i_forward_raw_hash"] = (
        out["forward_raw_center_hash_r1"]
    )

    out["i_forward_edit_hash"] = (
        out["forward_edit_signature_hash"]
    )

    out["i_forward_component_hash"] = (
        out["forward_component_evidence_hash"]
    )

    if store_readable:
        for src, dst in (
            (
                "raw_center_signature_r0",
                "forward_raw_center_signature_r0",
            ),
            (
                "raw_center_signature_r1",
                "forward_raw_center_signature_r1",
            ),
            (
                "raw_center_signature_r2",
                "forward_raw_center_signature_r2",
            ),
            (
                "edit_signature",
                "forward_edit_signature",
            ),
            (
                "component_evidence_signature",
                "forward_component_evidence_signature",
            ),
            (
                "component_topology_signature",
                "forward_component_topology_signature",
            ),
        ):
            if src in f:
                out[dst] = f[src]

    return out
