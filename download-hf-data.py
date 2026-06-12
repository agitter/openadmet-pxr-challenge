import pandas as pd

files = {
    "train": "pxr-challenge_TRAIN.csv",
    "test": "pxr-challenge_TEST_BLINDED.csv",
    "train_counter": "pxr-challenge_counter-assay_TRAIN.csv",
    "test_structure": "pxr-challenge_structure_TEST_BLINDED.csv",
    "train_single": "pxr-challenge_single_concentration_TRAIN.csv",
    "train_crudes": "pxr-challenge_htchem-libraries_TRAIN.csv",
    "train_semi_pure": "pxr-challenge_96-compound-uscale-semi-pure_TRAIN.csv",
    "test_phase1": "pxr-challenge_TEST_PHASE_1_UNBLINDED.csv",
}

for name, fname in files.items():
    df = pd.read_csv(f"hf://datasets/openadmet/pxr-challenge-train-test/{fname}")
    df.to_csv(f"data/{fname}", index=False)
    print(name, df.shape)
