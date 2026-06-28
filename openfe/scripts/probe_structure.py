"""Probe the full openfe result JSON structure to find convergence fields."""
import json, glob, os

# Find the -15 kcal/mol edge complex result
for d in glob.glob('openfe/production/146/rbfe_*OADMET-0006465*OCNT-2311117*complex*'):
    for sub in ['result.json', 'quickrun_output/result.json']:
        p = os.path.join(d, sub)
        if os.path.exists(p) and os.path.getsize(p) > 100:
            data = json.load(open(p))

            # 1. The protocol_result.data hashed entry
            pr_data = data.get('protocol_result', {}).get('data', {})
            for hashkey, val in pr_data.items():
                print('=== protocol_result.data entry ===')
                print('type:', type(val))
                if isinstance(val, list):
                    print('list len:', len(val))
                    if val and isinstance(val[0], dict):
                        print('first item keys:', list(val[0].keys()))
                elif isinstance(val, dict):
                    print('keys:', list(val.keys()))
                break

            # 2. All unit_result names and their output keys
            print()
            print('=== all unit_results ===')
            for k, v in data.get('unit_results', {}).items():
                name = v.get('name', '')
                outputs = v.get('outputs', {})
                okeys = list(outputs.keys()) if isinstance(outputs, dict) else str(type(outputs))
                print(f'  {name}')
                print(f'    output keys: {okeys}')
            break
