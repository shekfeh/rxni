"""Central semantic validation boundary.

Runtime v6.6 Stage E validation boundary.

All 82 registered directional FGI events have executable validators.  The
canonical Stage-E pipeline consumes ACCEPTED results only; REJECTED and
AMBIGUOUS hypotheses remain diagnostics, and NOT_RUN is a canonicalization
error.
"""

from __future__ import annotations

from typing import Any

from ..contracts import ContractRegistry, load_default_registry
from ..types import SemanticRecord, SemanticValidation, ValidationStatus
from ..version import FGI_CONTRACT_VERSION
from .priority1 import validate_priority1_directional
from .vc1_center_state import validate_vc1_center_state
from .wave_b_functional_centers import validate_wave_b_functional_centers
from .wave_c_bond_attachment import validate_wave_c_bond_attachment
from .wave_d_multicenter import validate_wave_d_multicenter


class ValidationEngine:
    """Validate detector hypotheses against explicit chemistry contracts.

    Unsupported registered events remain NOT_RUN. Absence of a native
    validator must never be interpreted as acceptance.
    """

    def __init__(self, registry: ContractRegistry | None = None):
        self.registry = registry or load_default_registry()

    def validate_directional(
        self,
        facts: dict[str, Any],
        record: SemanticRecord,
    ) -> SemanticValidation:
        contract = self.registry.get(record.event)
        if contract is None:
            return SemanticValidation(
                label=record.event,
                status=ValidationStatus.NOT_RUN,
                veto_reasons=("no_registered_contract",),
                contract_version=FGI_CONTRACT_VERSION,
            )

        result = validate_priority1_directional(facts, record)
        if result is not None:
            return result

        result = validate_vc1_center_state(
            facts,
            record,
            contract,
        )
        if result is not None:
            return result

        result = validate_wave_b_functional_centers(
            facts,
            record,
            contract,
        )
        if result is not None:
            return result

        result = validate_wave_c_bond_attachment(
            facts,
            record,
            contract,
        )
        if result is not None:
            return result

        result = validate_wave_d_multicenter(
            facts,
            record,
            contract,
        )
        if result is not None:
            return result

        return SemanticValidation(
            label=record.event,
            status=ValidationStatus.NOT_RUN,
            validator_class=contract.validator_class,
            veto_reasons=("native_validator_not_ported_yet",),
            contract_version=FGI_CONTRACT_VERSION,
        )
