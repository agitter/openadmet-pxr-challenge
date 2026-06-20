import os
from pathlib import Path
results = Path('results')
rows = []
for d in sorted(results.iterdir()):
    if not d.is_dir(): continue
    cid = d.name.replace('network_setup_', '')
    jsons = list((d / 'transformations').glob('*.json')) if (d / 'transformations').exists() else []
    has_network = (d / 'network_setup.json').exists()
    rows.append((cid, len(jsons), has_network))
    if not has_network or len(jsons) == 0:
        print(f'PROBLEM: cluster {cid}: {len(jsons)} JSONs, network_setup.json={has_network}')

total = sum(r[1] for r in rows)
print(f'\nTotal transformation JSONs: {total}')
print(f'Expected ~{len(rows)*2*2} (rough estimate: avg 4 edges x 2 legs x {len(rows)} clusters)')
# Distribution
counts = [r[1] for r in rows]
print(f'JSONs per cluster: min={min(counts)}, max={max(counts)}, mean={sum(counts)/len(counts):.1f}')
