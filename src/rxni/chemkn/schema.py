from __future__ import annotations
from dataclasses import dataclass, asdict
from typing import Any

ALIASES={
 "reaction_key":("reaction_key","query_reaction_key"),
 "template_hash":("retro_template_hash","template_hash"),
 "reaction_id":("reaction_id","member_id","id"),
 "catalyst":("catalyst","catalysts","catalyst_family"),
 "base":("base","bases","base_family"),
 "solvent":("solvent","solvents","solvent_family"),
 "reagent":("reagent","reagents","reagent_family"),
 "temperature":("temperature","temperature_c","temp_c","temp"),
 "time":("time","time_h","duration_h"),
 "condition_cluster":("condition_cluster","condition_card","condition_label"),
}

@dataclass
class ConditionRecord:
    reaction_key: str|None=None
    template_hash: str|None=None
    reaction_id: str|None=None
    catalyst: Any=None
    base: Any=None
    solvent: Any=None
    reagent: Any=None
    temperature: float|None=None
    time: float|None=None
    condition_cluster: str|None=None
    raw: dict[str,Any]|None=None
    def as_dict(self): return asdict(self)
