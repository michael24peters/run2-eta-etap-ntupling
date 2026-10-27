# README

To run, make sure you are in the base directory and type: `lb-run DaVinci/v45r8 ipython src/ana.py`. There are boolean values in `ana.py` that should be set for the particular circumstance (see Usage).

The ntuple is broken into the `tag`, `prt`, and `mc`. `tag` and `prt` contain reconstruction-level event information; `mc` contains generator-level event information. There are a series of indices which match reco-level events to gen-level events for the purposes of background analysis.

For local runs, a `.root` file will be generated in the `ntuples/` directory.

`src/` contains the files used to fill the ntuple and other intermediate operations.

## Usage

In `ana.py` you will see near the top the following flags:

- `IS_MC`: True = MC, False = Run 2 data
- `IS_SIGNAL`: True = signal, False = minbias
- `IS_SAMPLE`: True = local sample (for testing), False = analysis production
- `DECAY`: string of decay type. Current options:
    - `'eta2mumu'`
    - `'eta2mumugamma'`
    - `'eta2mumumumu'`
    - `'eta2mumuee'`

Change these flags as appropriate.
