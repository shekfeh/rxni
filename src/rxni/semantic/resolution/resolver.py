"""Resolution boundary for validated semantic claims.

The native resolver may:
- select primary ownership,
- suppress broader overlapping labels,
- retain compatible co-events,
- preserve multiplicity.

It MUST NOT invent chemistry, repair failed hypotheses, or derive redox.
The current v6.5.6 resolver remains reference-only until ownership rules are
ported with 509 parity tests.
"""
