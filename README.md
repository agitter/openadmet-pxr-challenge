# OpenADMET Predicting PXR Induction Blind Challenge
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21084637.svg)](https://doi.org/10.5281/zenodo.21084637)

An OpenFE-based approach to the [OpenADMET Predicting PXR Induction Blind Challenge](https://huggingface.co/spaces/openadmet/pxr-challenge).
The [writeup](writeup.md) describes the full methodology and results.

- `/claude`: Scripts and files downloaded from a Claude session
- `/data`: Data files from [Hugging Face](https://huggingface.co/datasets/openadmet/pxr-challenge-train-test)
- `/docking`: Files for docking with [GNINA](https://github.com/gnina/gnina)
- `/external`: Refined PXR structures as a [submodule](https://github.com/OpenADMET/pxr_xtal_re-refinement) and files from the [organizers](https://github.com/OpenADMET/PXR-Challenge-Tutorial/tree/main/evaluation).
- `/openfe`: Files for running [OpenFE](https://github.com/OpenFreeEnergy/openfe)
- `/submissions`: Scraped from [Hugging Face](https://openadmet-pxr-challenge.hf.space/config)

Most code was written or drafted by Claude Sonnet 4.6, Claude Opus 4.6, Claude Opus 4.8, and GPT-5.5 Instant.

## Citation
```
@article{gitter_openadmet_pxr_2026,
	title = {Structure-based drug discovery for the {OpenADMET} {Predicting} {PXR} {Induction} {Blind} {Challenge}},
	url = {https://github.com/agitter/openadmet-pxr-challenge},
	doi = {10.5281/zenodo.21084637},
	journal = {Zenodo},
	author = {Gitter, Anthony},
	month = jun,
	year = {2026},
}
```

## Docking notes
Initially only cluster representatives were run with GNINA.
To compare docking results with OpenFE results, a second round of docking tested all test compounds and selected training compounds (anchors), including reruning the previous cluster representatives.
The original run used numeric cluster ID subdirectories.
The second run used zero-padded ligand IDs prefixed with T (test) or A (training anchor).

## Preparing OpenFE
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
