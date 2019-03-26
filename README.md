# FlyTracker_Fixes
Fixes for computing relative features for more than 2 flies with FlyTracker.

DISCLAIMER: The code present in this repo represents a modification of the original code created by [Eyrún Eyjólfsdóttir](http://www.vision.caltech.edu/~eeyjolfs/) of the Computational Vision Lab at Caltech. Source code, along with other relevant and useful information can be gotten [here](http://www.vision.caltech.edu/Tools/FlyTracker/).

ISSUE: The current version of FlyTracker is able to compute relative features (distance between flies, angle between flies, facing angle, and leg distance) for only two flies. If more than two flies are present, FlyTracker will not compute these features nor, consequently, export the respective .m files into the JAABA folder. This is unfortunate, as these features provide important information with which machine learning algorithms such as [JAABA](http://jaaba.sourceforge.net/) can accurately predict and automatically annotate specific behaviors.

CONTENT:
- FlyTracker Patch

For detailed description of these tools and how to use them, head over the [Wiki section](https://github.com/miguelgaspar24/FlyTracker_Fixes/wiki)!
