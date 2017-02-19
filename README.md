# FlyTracker
Fixes for computing relative features for more than 2 flies.

DISCLAIMER: This code is a modification of the original source code created by Eyrún Eyjólfsdóttir of the Computational Vision Lab at Caltech (profile at: http://www.vision.caltech.edu/~eeyjolfs/). Source code, along with other relevant and useful information can be gotten at: http://www.vision.caltech.edu/Tools/FlyTracker/.

ISSUE: The current version of FlyTracker (see link for source code above) is able to compute relative features (distance between flies, angle between flies, facing angle, and leg distance) for only a single pair of flies. If more than two flies are present, FlyTracker will not compute these features them as .m files into the JAABA folder. This is unfortunate, as these features provide important information with which machine learning algorithms such as JAABA (Janelia Automatic Animal Behavior Annotator) can accurately predict and automatically annotate specific behaviors.

CONTENT:
- Fly Tracker MultiFix
- FlyTracker Patch


## FlyTracker MultiFix

Here you will find a single file, which is a reworked version of the 'feat_compute.m' file of the original FlyTracker. You should copy this file into your current FlyTracker->tracking folder, and replace the file of the same name that is already there. From that point on, anytime FlyTracker will run, it will automatically compute all relative features, irrespective of the number of flies present in the current region of interest. 

---------------------

## FlyTracker Patch

This patch is for those who have already run data through the original FlyTracker and are thus missing the relative features. If this is your case, then you DO NOT need to rerun your data through the new, MultiFix version described above. This patch will compute all relative features, generate .xls/.csv files that incorporate those newly computed relative features, and export all relative features and their normalizations as .m files ito the JAABA. In essence, it allows users to get the missing relative features if they have run the standard version of FlyTracker without having to rerun the whole tracker on all the data again. Because the patch functions only on the data that already exists from the initial tracker run, the process is much faster than rerunning the tracker on the data again, and the end result is the exact same as running FlyTracker MultiFix.

##### Mode of Use

- Open main_patch.m;
- Make sure to set 'folderspath' to the fold that contains your videos and FlyTracker output folders;
- You can now run the script, which will fetch the remaing scripts.

##### Notes

The 'helper' folder contains two debugging tools. 'Truncator' is used to remove the computed relatives features, in case one wants the file to revert back to its pre-patch form. This only affects the '*-feat.mat' file.
The 'matrix_comparison' is used to compare the personal and environmental features before and after the patch, in order to verify that patching did not alter previously existing data.
