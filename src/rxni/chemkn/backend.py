from __future__ import annotations

from typing import Any

from .evidence import ConditionEvidenceStore
from .ranking import summarize_conditions


class ChemKNBackend:
    """Evidence-first condition recommendation backend.

    Structural/semantic classification is upstream.  Runtime-v6.5.3 may pass a
    query-conditioned member order, but this backend never changes the reaction
    interpretation itself.
    """

    def __init__(self, store: ConditionEvidenceStore, *, min_records: int = 5):
        self.store = store
        self.min_records = int(min_records)

    @classmethod
    def from_path(cls, path, **kw):
        return cls(ConditionEvidenceStore.load(path), **kw)

    def recommend(self, retrieval: dict[str, Any], *, top_templates: int = 20) -> dict[str, Any]:
        templates = retrieval.get("templates") or []
        members = retrieval.get("members") or []
        qcond = retrieval.get("query_condition_evidence") or {}
        conditioned = bool(
            qcond.get("reagents") or qcond.get("solvents") or qcond.get("catalysts")
            or qcond.get("bases") or qcond.get("additives") or qcond.get("acids")
            or qcond.get("temperature_c") is not None or qcond.get("reaction_time_h") is not None
        )
        tiers = []; selected = []; fallback = None

        if templates:
            th = templates[0].get("template_hash")
            recs = self.store.for_template(th)
            tiers.append({"tier": 0, "name": "top_template", "template_hash": th, "n_records": len(recs)})
            if len(recs) >= self.min_records:
                selected = recs; fallback = 0

        if not selected:
            keys = [m["reaction_key"] for m in members if m.get("exact_edit")]
            recs = self.store.for_reactions(keys)
            tiers.append({"tier": 1, "name": "exact_edit_members", "n_members": len(keys), "n_records": len(recs)})
            if len(recs) >= self.min_records:
                selected = recs; fallback = 1

        if not selected:
            keys = [m["reaction_key"] for m in members[:max(50, top_templates * 5)]]
            recs = self.store.for_reactions(keys)
            tier_name = "ranked_query_conditioned_neighbors" if conditioned else "ranked_structural_neighbors"
            tiers.append({"tier": 2, "name": tier_name, "n_members": len(keys), "n_records": len(recs)})
            selected = recs; fallback = 2

        summary = summarize_conditions(selected)
        return {
            "schema_version": "condition_recommendation_v653",
            "fallback_level": fallback,
            "evidence_tiers": tiers,
            "query_condition_evidence": qcond,
            "ranking_mode": retrieval.get("ranking_mode"),
            "summary": summary,
            "evidence_records": [r.as_dict() for r in selected[:200]],
        }
