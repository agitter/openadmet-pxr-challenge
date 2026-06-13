# OpenADMET PXR Blind Challenge
- `/claude`: Scripts and files downloaded from a Claude session
- `/data`: Data files from HuggingFace
- `/docking`: Files for docking with GNINA
- `/external`: Refined PXR structures
- `/openfe`: Files for running OpenFE

Most code written or drafted by Claude Sonnet 4.6 and GPT-5.5 Instant.

The Apptainer command was run within an interactive Docker session using `ghcr.io/apptainer/apptainer:1.4.5`:
```
apptainer pull openfe_1.11.1.sif \
  oras://ghcr.io/openfreeenergy/openfe:1.11.1-apptainer
```
Upload to CHTC:
```
scp openfe_1.11.1.sif agitter@ap2001.chtc.wisc.edu:/staging/a/agitter/containers/
```
Test the image in an interactive HTCondor session:
```
apptainer run --nv /staging/a/agitter/containers/openfe_1.11.1.sif python -c "import sys; print(sys.version)"
```

```commandline
$ python analyze_pxr_ensemble.py
Found 64 structure directories
  1m13: ligand=HYF (39 heavy atoms, formula=C35H68O4), pocket_residues=67, resolution=None
  1nrl: ligand=SRL (33 heavy atoms, formula=C24H54O7P2), pocket_residues=62, resolution=None
  1skx: ligand=RFP (51 heavy atoms, formula=C38H71NO12), pocket_residues=64, resolution=None
  2o9i: ligand=444 (31 heavy atoms, formula=C17H26F9NO3S), pocket_residues=61, resolution=2.8
  2qnv: ligand=CDZ (29 heavy atoms, formula=C25H50O4), pocket_residues=49, resolution=None
  3hvl: ligand=SRL (33 heavy atoms, formula=C24H54O7P2), pocket_residues=64, resolution=None
  3r8d: ligand=PNU (20 heavy atoms, formula=C13H25ClN4OS), pocket_residues=47, resolution=None
  4ny9: ligand=2Q4 (30 heavy atoms, formula=C23H45ClN2O4), pocket_residues=61, resolution=None
  4s0t: ligand=40U (29 heavy atoms, formula=C23H43ClN2O3), pocket_residues=61, resolution=None
  4x1f: ligand=3WF (22 heavy atoms, formula=C20H34O2), pocket_residues=50, resolution=None
  4x1g: ligand=3WF (22 heavy atoms, formula=C20H34O2), pocket_residues=48, resolution=None
  4xhd: ligand=40U (29 heavy atoms, formula=C23H43ClN2O3), pocket_residues=60, resolution=None
  5a86: ligand=D7E (28 heavy atoms, formula=C18H33ClF3N3O2S), pocket_residues=50, resolution=None
  5x0r: ligand=4WH (30 heavy atoms, formula=C22H45N3O4S), pocket_residues=64, resolution=None
  6bns: ligand=XGH (39 heavy atoms, formula=C23H39F7N2O5S2), pocket_residues=73, resolution=None
  6dup: ligand=HCJ (28 heavy atoms, formula=C20H33F3N2O3), pocket_residues=53, resolution=None
  6hj2: ligand=P06 (35 heavy atoms, formula=C23H44F3N5O2S2), pocket_residues=65, resolution=None
  6hty: ligand=GRH (29 heavy atoms, formula=C21H41ClN2O4S), pocket_residues=61, resolution=2.22
  6nx1: ligand=L7D (31 heavy atoms, formula=C20H31F7O3S), pocket_residues=53, resolution=2.27
  6p2b: ligand=NQD (31 heavy atoms, formula=C27H52O4), pocket_residues=60, resolution=None
  6s41: ligand=KUB (28 heavy atoms, formula=C17H31ClF3N3O2S2), pocket_residues=60, resolution=2.7
  6tfi: ligand=N6H (30 heavy atoms, formula=C22H40ClN3O4), pocket_residues=58, resolution=None
  6xp9: ligand=QCG (37 heavy atoms, formula=C29H52FN3O4), pocket_residues=57, resolution=2.27
  7ax9: ligand=S6H (18 heavy atoms, formula=C10H8Cl8), pocket_residues=49, resolution=None
  7axa: ligand=CL6 (25 heavy atoms, formula=C22H39ClN2), pocket_residues=53, resolution=None
  7axb: ligand=S68 (19 heavy atoms, formula=C9H10Cl6O3S), pocket_residues=46, resolution=None
  7axc: ligand=27H (26 heavy atoms, formula=C22H40O4), pocket_residues=48, resolution=None
  7axd: ligand=S6W (26 heavy atoms, formula=C12H20Cl2F6N4OS), pocket_residues=55, resolution=None
  7axe: ligand=S6Z (22 heavy atoms, formula=C15H28Cl2N2O3), pocket_residues=55, resolution=None
  7axf: ligand=S6T (21 heavy atoms, formula=C17H34ClNO2), pocket_residues=48, resolution=None
  7axg: ligand=TBY (13 heavy atoms, formula=C12H28Sn), pocket_residues=42, resolution=2.7
  7axh: ligand=27J (23 heavy atoms, formula=C18H34O5), pocket_residues=49, resolution=None
  7axi: ligand=EST (20 heavy atoms, formula=C18H30O2), pocket_residues=47, resolution=2.15
  7axj: ligand=CL6 (25 heavy atoms, formula=C22H39ClN2), pocket_residues=46, resolution=2.3
  7axk: ligand=EST (20 heavy atoms, formula=C18H30O2), pocket_residues=45, resolution=2.0
  7axl: ligand=EST (20 heavy atoms, formula=C18H30O2), pocket_residues=45, resolution=2.5
  7n2a: ligand=07F (31 heavy atoms, formula=C26H47FN2O2), pocket_residues=58, resolution=None
  7rio: ligand=5YX (40 heavy atoms, formula=C26H51F3N4O4S3), pocket_residues=70, resolution=None
  7riu: ligand=5YU (32 heavy atoms, formula=C21H42FN5O2S3), pocket_residues=59, resolution=2.05
  7riv: ligand=P06 (35 heavy atoms, formula=C23H44F3N5O2S2), pocket_residues=65, resolution=None
  7yfk: ligand=G3L (37 heavy atoms, formula=C28H46O9), pocket_residues=60, resolution=2.1
  8cct: ligand=UAI (19 heavy atoms, formula=C15H26Cl2O2), pocket_residues=48, resolution=None
  8cf9: ligand=UK6 (22 heavy atoms, formula=C20H38O2), pocket_residues=48, resolution=None
  8ch8: ligand=ULC (23 heavy atoms, formula=C18H34N2O2S), pocket_residues=54, resolution=2.15
  8e3n: ligand=VA0 (50 heavy atoms, formula=C37H69NO12), pocket_residues=69, resolution=None
  8eqz: ligand=WQB (32 heavy atoms, formula=C21H37F6NO3S), pocket_residues=58, resolution=None
  8f5y: ligand=JQ1 (31 heavy atoms, formula=C23H43ClN4O2S), pocket_residues=64, resolution=2.15
  8fpe: ligand=Y5B (39 heavy atoms, formula=C28H47F6NO3S), pocket_residues=58, resolution=2.3
  8r00: ligand=XFQ (17 heavy atoms, formula=C15H28O2), pocket_residues=43, resolution=1.95
  8r81: ligand=Y8B (29 heavy atoms, formula=C23H45N3O2S), pocket_residues=63, resolution=None
  8r82: ligand=Y7Q (37 heavy atoms, formula=C30H60N4O2S), pocket_residues=68, resolution=None
  8svo: ligand=WSO (34 heavy atoms, formula=C27H54N4O3), pocket_residues=70, resolution=None
  8svp: ligand=WSX (36 heavy atoms, formula=C28H56N4O4), pocket_residues=72, resolution=None
  8svq: ligand=WSX (36 heavy atoms, formula=C28H56N4O4), pocket_residues=70, resolution=None
  8svr: ligand=WT1 (38 heavy atoms, formula=C29H58N4O5), pocket_residues=74, resolution=None
  8svs: ligand=WU2 (37 heavy atoms, formula=C29H58N4O4), pocket_residues=71, resolution=None
  8svt: ligand=WU6 (36 heavy atoms, formula=C26H52N4O6), pocket_residues=72, resolution=None
  8svx: ligand=WU6 (36 heavy atoms, formula=C26H52N4O6), pocket_residues=72, resolution=None
  8szv: ligand=X1D (40 heavy atoms, formula=C28H47F6NO4S), pocket_residues=67, resolution=2.2
  9beq: ligand=AP1 (50 heavy atoms, formula=C38H71NO11), pocket_residues=69, resolution=2.6
  9fzg: ligand=A1IHA (60 heavy atoms, formula=None), pocket_residues=69, resolution=2.0
  9fzh: ligand=A1IHB (50 heavy atoms, formula=None), pocket_residues=66, resolution=2.5
  9fzi: ligand=SRL (33 heavy atoms, formula=C24H54O7P2), pocket_residues=63, resolution=None
  9fzj: ligand=SRL (33 heavy atoms, formula=C24H54O7P2), pocket_residues=62, resolution=None

Wrote pxr_structure_inventory.csv (64 rows)
Wrote pxr_pocket_residues.json

Structures with successfully extracted ligand SMILES: 62 / 64

NOTE: SMILES extracted from PDB connectivity via RDKit can have WRONG BOND ORDERS / protonation states (PDB format dt reliably encode these). Treat 'ligand_smiles' as a rough identity hint for similarity matching, not a validated ste. Cross-check important ligands manually against PDB ligand IDs (ligand_resname) at https://www.rcsb.org/ligand/<re
(base)
```