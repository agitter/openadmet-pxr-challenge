# Structure-based drug discovery for the OpenADMET Predicting PXR Induction Blind Challenge
Anthony Gitter  
Archived at TBD

## Introduction
The [OpenADMET Predicting PXR Induction Blind Challenge](https://huggingface.co/spaces/openadmet/pxr-challenge) provided a rich dataset that reflects multiple stages of hit discovery and lead optimization: a primary screen, dose-response screen, counter-screen, and dose-response screen on an analog expansion set.
The goal was to use any of that data to predict pEC50 (-log of EC50) for the 513 compound analog expansion set that served as the test set.
My strategy was formed based on a few goals, assumptions and constraints.
I started late on June 12, just two and a half weeks before the deadline.
This was after phase 1 ended, so I did not have opportunities for feedback on the live leaderboard, but I did have the unblinded analog set 1 data to use.
Because the target compounds were intentionally structurally similar to the 63 active compounds, I was skeptical of how well supervised learning approaches could model the precise quantitative structure activity relationships.
Finally, I was eager to use the challenge to try something different that is outside my research group's typical workflow.
I did this with [TabPFN](https://github.com/agitter/asap-polaris-admet-challenge/blob/main/writeup.md) on the [ASAP Discovery x OpenADMET Challenge](https://doi.org/10.1021/acs.jcim.5c02106) and [AI scientists](https://github.com/agitter/openadmet-expansionrx-challenge/blob/main/writeup.md) on the [OpenADMET + ExpansionRx Blind Challenge](https://huggingface.co/spaces/openadmet/OpenADMET-ExpansionRx-Challenge).

That combination of factors motivated me to select [OpenFE](https://docs.openfree.energy/en/stable/index.html) relative binding free energy (RBFE) simulation as my core methodology.
[Alyssa Travitz](https://alyssatravitz.com/) from the [OpenFE team](https://openfree.energy/team/) recently gave a seminar here at the University of Wisconsin-Madison that hooked me.
To preview my results, the rest of this report could be characterized as what my clinical biostatistics colleagues call a [_futility analysis_](https://pubmed.ncbi.nlm.nih.gov/17128426/).
It became clear pre-submission that my approach had failed.
However, the process was educational.
I confirmed OpenFE is a great match for the computing infrastructure I have available, and I could see myself returning to RBFE analysis again (in closer consultation with my molecular simulation colleagues).

## Methods
Outline to be expanded:
- Select refined PXR structures
- Run [GNINA](https://github.com/gnina/gnina) docking
- Run OpenFE RBFE simulations
- Analyze docking and RBFE results and RBFE failures
- Calibrate docking and RBFE scores and sweep over ensemble models
- Select a final ensemble model and score test compounds

Most code was written or drafted by Claude Sonnet 4.6, Claude Opus 4.6, or Claude Opus 4.8.
I also used GPT-5.5 Instant to prepare and test the OpenFE Singularity image and scripts.

## Results
Text to come

![Docking scores versus pEC50](docking/docking_analysis_extended/docking_vs_pec50_correlation.png)

![Comparisons of docking score, RBFE score, and pEC50](openfe/method_comparison.png)

![Assessment of ensemble model options](openfe/model_tuning_diagnostics.png)

![Threshold selection for simple docking and RBFE ensemble](openfe/blend_weight_selection.png)

![Submission score inspections](openfe/submission_visualizations.png)

![Summary of GPU computing for RBFE simulations](openfe/compute_accounting.png)

![GPU computing timeline](openfe/compute_timeline.png)

## Discussion
- Summarize the auto-generated limitations document