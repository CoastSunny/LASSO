# Sparse EEG/MEG source estimation via a group lasso
Michael Lim, Justin M. Ales, Benoit R. Cottereau, Trevor Hastie, Anthony M. Norcia,

This repository contains the code to recreate Figure 7,9,10 from the paper.


Running the code required downloading the associated data from: https://osf.io/4mz47/

Extract the data into a directory named "datafiles" containing the subdirectories: "anatomy", "eegdata", "forwardSolutions".  If the datafiles directory is correctly placed in the same folder as mrLASSO.m it will be found automatically. If not the code will ask you to select the directory.

From MATLAB run: mrLASSO

The code will ask several questions.  First whether you want to compute the minimum norm solution or use the precomputed solution provided. The cross participant minimum-norm solution takes a lengthy amount of time (~20 minutes at last test). Therefore, we've provided a precumputed solution along with the data files.   After generating several data figures the next thing that will be asked is whether or not to plot topographies. If you choose to plot topographies you will be asked if you want to plot a single participant or all participants.  As each topography is plotted on a high-resolution scalp mesh they can be intensive to render.  Plotting a single participant will allow interacting in 3d. Plotting all participants will render each participant, take a snapshot, then create a view of all participants together as in Figure 10.


