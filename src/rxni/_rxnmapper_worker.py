"""Standalone RXNMapper worker; executable using an optional Python environment.

No RXNI imports: this worker can run from an independent interpreter with
`rxnmapper` installed. The parent application validates the returned mapping.
"""
import json
import sys


def main():
    from rxnmapper import RXNMapper  # intentionally optional / isolated
    data = json.load(sys.stdin)
    mapper = RXNMapper()
    result = mapper.get_attention_guided_atom_maps([data['reaction']])[0]
    print(json.dumps({'mapped_reaction':result['mapped_rxn'],
                      'confidence':result.get('confidence')}))


if __name__ == '__main__':
    main()
