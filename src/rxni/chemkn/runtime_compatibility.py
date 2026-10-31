from __future__ import annotations
from .knowledge_schema import RuntimeCompatibilityQuery


def from_retrieval_result(result: dict) -> RuntimeCompatibilityQuery:
    record=result.get("query_record") or {}
    sem=record.get("semantics") or {}
    events=set()
    for bucket in ("primary_events","semantic_events","specific_events","local_events","primitives","topology_evidence","directional_conversions"):
        for x in sem.get(bucket) or []:
            if not isinstance(x,dict): continue
            for k in ("event","subtype","directional_conversion"):
                if x.get(k): events.add(str(x[k]))
    primary=result.get("primary_event")
    if primary: events.add(str(primary))
    formed=sem.get("formed_bond_classes") or []
    for b in formed:
        key=str(b).upper().replace("-","_")
        if key=="C_O": events.add("c_o_formation")
        elif key=="C_N": events.add("c_n_formation")
        elif key=="C_C": events.add("c_c_formation")
        elif key=="C_S": events.add("c_s_formation")
    topo=sem.get("event_topology") or {}
    for b in topo.get("formed_bond_classes") or []:
        key=str(b).upper().replace("-","_")
        if key=="C_O": events.add("c_o_formation")
        elif key=="C_N": events.add("c_n_formation")
        elif key=="C_C": events.add("c_c_formation")
        elif key=="C_S": events.add("c_s_formation")
    # Direction stays unset unless runtime exposes the same controlled vocabulary
    # as the compatibility overlays; guessing here would make hard gates unsafe.
    return RuntimeCompatibilityQuery.from_values(events=events)
