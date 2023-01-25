
## write workspace
python SIGfits.py (set rootFit and functional form, fill the workspace)

SIGfits.py call the supporting files are:
* prepareFits.py (changed directory, working MVA points, preselection)
* LoadTree.py

## write DATACARDS and compute limits
./combineCommand.sh

## make limit Plots
python limitPlot.py