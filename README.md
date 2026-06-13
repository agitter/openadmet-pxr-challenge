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
