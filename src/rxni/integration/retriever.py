"""RXNI reaction retriever: stable public alias, frozen retrieval algorithm.

`ReactionRetriever` aliases the historical `V66Retriever` class. Older class
and scoring names are compatibility/provenance identifiers, not statements
about the scientific authority of any chosen library. Query records are built
by the installed RXNI semantic engine; stored profiles are read from the
configured V2 representation provider and can have their own independently
attested generation lineage. Preserve all stored keys/hashes/scoring contracts.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
from typing import Any
import hashlib
import sqlite3
import time

# Query record construction is the semantic routing seam: use the current
# Runtime-v6.7 semantic engine while leaving the retrieval/scoring algorithms
# frozen at v6.6.
from rxni.retrieval.record import build_retrieval_record

from rxni.retrieval import (
    member_scores,
)
from rxni.retrieval.storage.provider import (
    V2StoredRecordProvider,
)
from rxni.method.reagents import (
    extract_reagents_from_reaction,
)
from rxni.integration.method_orchestrator import (
    orchestrate_post_retrieval,
    _semantic_primary_event,
    _score_value,
)

from rxni.retrieval import derive_query_structural_keys
from rxni.retrieval.libraries import resolve_library
from rxni.retrieval.libraries.foundation_schema import (
    load_foundation_profile, quote_ident,
)

CHANNEL_TO_COLUMN = {'exact_edit': 'forward_edit_signature_hash', 'raw_r0': 'forward_raw_center_hash_r0', 'raw_r1': 'forward_raw_center_hash_r1', 'raw_r2': 'forward_raw_center_hash_r2', 'component_topology': 'forward_component_topology_hash', 'component_evidence': 'forward_component_evidence_hash'}
SPECIFICITY_PRIORITY = ('raw_r2', 'raw_r1', 'exact_edit', 'component_evidence', 'raw_r0', 'component_topology')


def _query_key(reaction: str) -> str:
    # Output/provenance identity only. Candidate identity is always the
    # library's stored reaction_key.
    return hashlib.sha256(str(reaction).encode("utf-8")).hexdigest()


def _sqlite_text(value: Any) -> str:
    # Structural-index values are persisted as exact text identity.
    # numpy / Arrow integer scalars stringify losslessly here.
    return str(value)


class _MetadataLookup:
    """Schema-discovered read-only reaction/member metadata lookup."""

    REQUIRED = {
        "reaction_key",
        "forward_reaction_smiles",
    }

    PREFERRED = (
        "reaction_id",
        "forward_reaction_smiles",
        "retro_template_hash",
        "retro_template_number",
        "membership_occurrence_count",
    )

    def __init__(self, path: str | Path, profile=None):
        self.path = Path(path)
        if not self.path.exists():
            raise FileNotFoundError(self.path)
        self.conn = sqlite3.connect(
            f"file:{self.path}?mode=ro",
            uri=True,
        )
        if profile is None:
            self.table, self.columns = self._discover_table()
            self.mapping = {name: name for name in self.columns}
        else:
            self.table = profile.metadata_table
            self.mapping = dict(profile.metadata_columns)
            self.columns = set(self.mapping)

    def _discover_table(self):
        tables = [
            str(x[0])
            for x in self.conn.execute(
                "SELECT name FROM sqlite_master "
                "WHERE type='table' ORDER BY name"
            ).fetchall()
        ]

        candidates = []
        for table in tables:
            cols = {
                str(x[1])
                for x in self.conn.execute(
                    f'PRAGMA table_info("{table}")'
                ).fetchall()
            }
            if not self.REQUIRED.issubset(cols):
                continue
            score = sum(1 for x in self.PREFERRED if x in cols)
            candidates.append((score, table, cols))

        if not candidates:
            raise RuntimeError(
                "No runtime_lookup table contains required fields: "
                + ", ".join(sorted(self.REQUIRED))
            )

        candidates.sort(key=lambda x: (-x[0], x[1]))
        _, table, cols = candidates[0]
        return table, cols

    def get_many(self, reaction_keys):
        keys = [str(x) for x in reaction_keys]
        if not keys:
            return {}

        cols = [
            "reaction_key",
            *[x for x in self.PREFERRED if x in self.columns],
        ]
        # reaction_key may also be in PREFERRED in future; dedupe.
        cols = list(dict.fromkeys(cols))

        out = {}
        chunk = 800
        for i in range(0, len(keys), chunk):
            block = keys[i:i + chunk]
            q = ",".join("?" for _ in block)
            sql = (
                "SELECT "
                + ",".join(quote_ident(self.mapping[x]) + ' AS ' + quote_ident(x) for x in cols)
                + ' FROM ' + quote_ident(self.table) + ' '
                + 'WHERE ' + quote_ident(self.mapping['reaction_key']) + f' IN ({q})'
            )
            for row in self.conn.execute(sql, block).fetchall():
                obj = dict(zip(cols, row))
                out[str(obj["reaction_key"])] = obj
        return out

    def close(self):
        try:
            self.conn.close()
        except Exception:
            pass


def _aggregate_templates(
    members: list[dict[str, Any]],
    *,
    mode: str,
) -> list[dict[str, Any]]:
    """Native dict-based equivalent of the frozen template aggregation rule."""
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for member in members:
        th = member.get("template_hash")
        if th:
            grouped[str(th)].append(member)

    rows = []
    for template_hash, group in grouped.items():
        ordered = sorted(
            group,
            key=lambda m: (
                -float(_score_value(m.get("scores") or {}, mode)),
                str(m.get("reaction_key") or ""),
            ),
        )
        vals = [
            float(_score_value(m.get("scores") or {}, mode))
            for m in ordered
        ]
        top3 = vals[:3]
        best = vals[0] if vals else 0.0
        mean_top3 = (
            sum(top3) / len(top3)
            if top3 else 0.0
        )
        score = 0.75 * best + 0.25 * mean_top3
        rows.append({
            "template_hash": template_hash,
            "template_number": ordered[0].get("template_number"),
            "score": score,
            "best_member_score": best,
            "mean_top3_member_score": mean_top3,
            "member_count": len(group),
            "best_member_key": ordered[0].get("reaction_key"),
            "best_member_reaction": ordered[0].get("reaction"),
        })

    rows.sort(
        key=lambda x: (
            -float(x["score"]),
            str(x["template_hash"]),
        )
    )
    return rows


class V66Retriever:
    """Frozen v6.6 retrieval machinery with a Runtime-v6.7 query semantic record."""

    def __init__(
        self,
        library: str | Path,
        *,
        row_group_cache_size: int = 4,
    ):
        # Storage layout is an adapter concern, not a chemical or scoring rule.
        # The old five-file directory ABI remains the default read path.
        resolved = resolve_library(library)
        self.library = resolved.root
        self.structural_db = resolved.structural_index
        self.lookup_db = resolved.metadata_lookup
        self.representations = resolved.representations
        self.locator = resolved.locator
        self.manifest = resolved.retrieval_manifest
        self.foundation_profile = (load_foundation_profile(resolved.foundation_profile)
                                   if resolved.foundation_profile is not None else None)
        self.structural_table = (self.foundation_profile.structural_table
                                 if self.foundation_profile else 'structural_index')
        self.structural_mapping = (self.foundation_profile.structural_columns
                                   if self.foundation_profile else {c:c for c in
                                       ('reaction_key', *CHANNEL_TO_COLUMN.values())})
        self.lookup = _MetadataLookup(self.lookup_db, profile=self.foundation_profile)
        self.provider = V2StoredRecordProvider(
            self.representations,
            self.locator,
            self.manifest,
            row_group_cache_size=row_group_cache_size,
        )

    def close(self):
        self.lookup.close()

    def structural_candidates(
        self,
        keys: dict[str, Any],
        *,
        per_channel_limit: int = 5000,
    ) -> tuple[dict[str, list[str]], dict[str, int]]:
        conn = sqlite3.connect(
            f"file:{self.structural_db}?mode=ro",
            uri=True,
        )
        out = {}
        available = {}
        try:
            for channel, col in CHANNEL_TO_COLUMN.items():
                val = keys.get(col)
                if val is None:
                    out[channel] = []
                    available[channel] = 0
                    continue

                sval = _sqlite_text(val)
                physical_col = quote_ident(self.structural_mapping[col])
                physical_key = quote_ident(self.structural_mapping['reaction_key'])
                table = quote_ident(self.structural_table)
                available[channel] = int(
                    conn.execute(
                        f'SELECT COUNT(*) FROM {table} '
                        f'WHERE {physical_col}=?',
                        (sval,),
                    ).fetchone()[0]
                )
                rows = conn.execute(
                    f'SELECT {physical_key} FROM {table} '
                    f'WHERE {physical_col}=? LIMIT ?',
                    (sval, int(per_channel_limit)),
                ).fetchall()
                out[channel] = [str(x[0]) for x in rows]
        finally:
            conn.close()

        return out, available

    @staticmethod
    def _select_candidates(
        channel_map: dict[str, list[str]],
        *,
        max_members: int,
        policy: str,
    ):
        channels_by_key = defaultdict(set)
        for channel, keys in channel_map.items():
            for key in keys:
                channels_by_key[key].add(channel)

        if policy == "legacy":
            selected = sorted(
                channels_by_key,
                key=lambda key: (
                    -len(channels_by_key[key]),
                    key,
                ),
            )[:max_members]
        elif policy == "specificity":
            selected = []
            seen = set()
            for channel in SPECIFICITY_PRIORITY:
                for key in channel_map.get(channel, []):
                    if key in seen:
                        continue
                    selected.append(key)
                    seen.add(key)
                    if len(selected) >= max_members:
                        break
                if len(selected) >= max_members:
                    break
        else:
            raise ValueError(
                f"Unknown candidate policy: {policy!r}"
            )

        selected_channel_counts = {
            channel: 0
            for channel in channel_map
        }
        for key in selected:
            for channel in channels_by_key[key]:
                selected_channel_counts[channel] += 1

        return selected, selected_channel_counts, channels_by_key

    def retrieve(
        self,
        reaction: str,
        *,
        mode: str = "sem_multi",
        top_k: int = 20,
        max_members: int = 3000,
        per_channel_limit: int = 5000,
        candidate_policy: str = "specificity",
        reagents: list[str] | None = None,
        solvents: list[str] | None = None,
        catalysts: list[str] | None = None,
        bases: list[str] | None = None,
        additives: list[str] | None = None,
        acids: list[str] | None = None,
        temperature_c: float | None = None,
        reaction_time_h: float | None = None,
        atmosphere: str | None = None,
        pressure: float | str | None = None,
        orderly_record: Any | None = None,
        reagent_weight: float = 0.20,
        profile_methods: bool = False,
        exclude_reaction_key: str | None = None,
        exclude_template: str | None = None,
        allowed_reaction_keys: set[str] | None = None,
        member_limit: int | None = None,
    ) -> dict[str, Any]:
        t0 = time.perf_counter()
        timings = {}
        qkey = _query_key(reaction)

        t = time.perf_counter()
        keys = derive_query_structural_keys(
            reaction,
            store_readable=False,
        )
        qrec = build_retrieval_record(
            reaction,
            record_id=qkey,
        )
        timings["query_representation"] = (
            time.perf_counter() - t
        )

        primary_event = _semantic_primary_event(qrec)

        t = time.perf_counter()
        channel_map, available_counts = (
            self.structural_candidates(
                keys,
                per_channel_limit=per_channel_limit,
            )
        )
        timings["candidate_generation"] = (
            time.perf_counter() - t
        )

        if exclude_reaction_key:
            for channel in channel_map:
                channel_map[channel] = [
                    key
                    for key in channel_map[channel]
                    if key != exclude_reaction_key
                ]

        if allowed_reaction_keys is not None:
            allowed = set(allowed_reaction_keys)
            for channel in channel_map:
                channel_map[channel] = [
                    key
                    for key in channel_map[channel]
                    if key in allowed
                ]

        (
            candidate_keys,
            selected_channel_counts,
            channels_by_key,
        ) = self._select_candidates(
            channel_map,
            max_members=max_members,
            policy=candidate_policy,
        )

        t = time.perf_counter()
        meta = self.lookup.get_many(candidate_keys)
        if exclude_template:
            candidate_keys = [
                key
                for key in candidate_keys
                if (
                    meta.get(key, {}).get(
                        "retro_template_hash"
                    )
                    != exclude_template
                )
            ]
        timings["metadata_lookup"] = (
            time.perf_counter() - t
        )

        t = time.perf_counter()
        self.provider.prefetch(candidate_keys)
        timings["representation_prefetch"] = (
            time.perf_counter() - t
        )

        t = time.perf_counter()
        base_members = []
        missing_provider = 0
        missing_metadata = 0

        for key in candidate_keys:
            m = meta.get(key)
            if not m:
                missing_metadata += 1
                continue

            crec = self.provider.get(
                record_id=key,
            )
            if crec is None:
                missing_provider += 1
                continue

            scores = dict(
                member_scores(qrec, crec)
            )
            candidate_primary = (
                _semantic_primary_event(crec)
            )

            reaction_smiles = str(
                m["forward_reaction_smiles"]
            )
            candidate_components = (
                extract_reagents_from_reaction(
                    reaction_smiles
                )
            )

            base_members.append({
                "reaction_key": key,
                "reaction": reaction_smiles,
                "template_hash": m.get(
                    "retro_template_hash"
                ),
                "template_number": m.get(
                    "retro_template_number"
                ),
                "membership_occurrence_count": (
                    m.get(
                        "membership_occurrence_count"
                    )
                    or 1
                ),
                "generation_channels": sorted(
                    channels_by_key.get(
                        key,
                        set(),
                    )
                ),
                "scores": scores,
                "exact_edit": (
                    "exact_edit"
                    in channels_by_key.get(
                        key,
                        set(),
                    )
                ),
                "metadata": {
                    "reaction_id": m.get(
                        "reaction_id"
                    ),
                    "candidate_primary_event":
                        candidate_primary,
                    "reagents":
                        candidate_components,
                },
            })

        timings["member_scoring"] = (
            time.perf_counter() - t
        )

        reagent_weight = min(
            1.0,
            max(0.0, float(reagent_weight)),
        )
        chemistry_weight = 1.0 - reagent_weight

        t = time.perf_counter()
        post = orchestrate_post_retrieval(
            base_members,
            reaction=reaction,
            primary_event=primary_event,
            mode=mode,
            chemistry_weight=chemistry_weight,
            reagent_weight=reagent_weight,
            cli_reagents=reagents,
            solvents=solvents,
            catalysts=catalysts,
            bases=bases,
            additives=additives,
            acids=acids,
            temperature_c=temperature_c,
            reaction_time_h=reaction_time_h,
            atmosphere=atmosphere,
            pressure=pressure,
            orderly_record=orderly_record,
        )
        timings["condition_method_orchestration"] = (
            time.perf_counter() - t
        )

        ranking_mode = post["ranking_mode"]
        members = list(post["members"])

        # M7E2 sorts condition-aware queries. For an unconditioned query,
        # preserve the historical explicit ranking step here.
        members.sort(
            key=lambda member: (
                -float(
                    _score_value(
                        member.get("scores") or {},
                        ranking_mode,
                    )
                ),
                str(
                    member.get(
                        "reaction_key"
                    )
                    or ""
                ),
            )
        )

        if not profile_methods:
            post["query_method_predictions"] = []
            post["global_method_landscape"] = None
            post["conditioned_method_landscape"] = None
            post["method_profile"] = None
            for member in members:
                md = member.get("metadata") or {}
                md.pop("method_predictions", None)
                md.pop("method_labels", None)

        candidate_category_counts = Counter(
            str(
                (m.get("metadata") or {}).get(
                    "candidate_category"
                )
                or "unclassified"
            )
            for m in members
        )

        templates = _aggregate_templates(
            members,
            mode=ranking_mode,
        )[:top_k]

        serialized_limit = (
            len(members)
            if member_limit == 0
            else (
                int(member_limit)
                if member_limit is not None
                else max(top_k * 5, 100)
            )
        )

        total = time.perf_counter() - t0
        timings["total"] = total

        return {
            "runtime_version": "6.7",
            "retrieval_facade":
                "runtime_v6_7.integration.retriever.V66Retriever",
            "semantic_runtime": "runtime_v6_7",
            "retrieval_generation": "runtime_v6_6_retrieval_v2",
            "query_reaction_key": qkey,
            "structural_keys": keys,
            "query_record": qrec,
            "primary_event": primary_event,
            "generation_counts": {
                key: len(value)
                for key, value in channel_map.items()
            },
            "generation_available_counts":
                available_counts,
            "candidate_policy":
                candidate_policy,
            "candidate_budget":
                max_members,
            "selected_channel_counts":
                selected_channel_counts,
            "selected_candidate_keys":
                len(candidate_keys),
            "candidate_member_count":
                len(members),
            "candidate_categories": [
                {
                    "category": key,
                    "count": value,
                    "fraction": (
                        value / len(members)
                        if members else 0.0
                    ),
                }
                for key, value
                in candidate_category_counts.most_common()
            ],
            "missing_provider":
                missing_provider,
            "missing_metadata":
                missing_metadata,
            "bulk_provider":
                self.provider.snapshot(),
            "query_condition_evidence":
                post["query_condition_evidence"],
            "query_reagents":
                post["query_reagents"],
            "query_method_predictions":
                post["query_method_predictions"],
            "reagent_weight":
                post["reagent_weight"],
            "condition_weight":
                post["condition_weight"],
            "ranking_mode":
                ranking_mode,
            "global_method_landscape":
                post["global_method_landscape"],
            "conditioned_method_landscape":
                post["conditioned_method_landscape"],
            "method_profile":
                post["method_profile"],
            "serialized_member_limit":
                serialized_limit,
            "members":
                members[:serialized_limit],
            "templates":
                templates,
            "mode":
                mode,
            "timings":
                timings,
            "elapsed_seconds":
                total,
        }


# Stable public class alias; implementation and response ABI unchanged.
ReactionRetriever = V66Retriever
