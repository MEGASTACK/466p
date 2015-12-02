# T10-Gait

Foot segment identification from footsteps on pressure pads

Google Drive folder is here: https://drive.google.com/open?id=0B2WhXU_IjE68cjlmMzZZdV9iMm8

Proposal is here: https://drive.google.com/open?id=1vuIcVz8KHbXqyjY4PWqkA_2YfWaEg1bcKuxf7oa9vcg

### Conditional Random Field training

File `main.m` trains a Conditional Random Field (CRF) model. It is heavily based off of `sophisticated_example.m`.

Setup:

1. Download JustinsGraphicalModelToolbox from the Google Drive, or from Slack. 
Unzip it and move it into the root of the project. The folder "JustinsGraphicalModelToolboxPublic" should be
in the same folder as this README.
2. Rename `JustinsGraphicalModelsToolboxPublic/External/toolbox/external/other/savefig.m` to
`JustinsGraphicalModelsToolboxPublic/External/toolbox/external/other/savefig_m.m`. We rename it
because it has the same name as a MATLAB builtin and causes problems.
3. Add the graphical model toolbox to MATLAB's path, recursively. Do this by
running `>> addpath(genpath('JustinsGraphicalModelsToolboxPublic'));` in MATLAB.

Train:

1. In MATLAB, run `>> main`.

### Semi-Automated Foot Labelling

File `labelling_automation/label_foot.m` is a function which will guide you 
through converting an image of a labelled foot into proper labelled foot data
which the system can consume in model training.

The paths are hardcoded. It must be run from within `labelling_automation`.

Example: labeling FinalData/NP40/3108

```
>> cd labelling_automation
>> label_foot('NP40', '3108'`)
```

In this case, the files `FinalData/NP40/3108.jpg` and `FinalData/NP40/3108.lst` must exist.
The jpg file is an image of an expert-segmented image, and the lst file is the footstep data
file for the given image.

The first step shows a picture of the labelled foot and instructs the user
to resize a rectangle to fit the outline of the colored squares. Each step after that
shows a picture of the labelled foot and instructs the user to trace one labelled region.

The region mask is then converted from pixel-coordinates to .lst-file coordinates and
the result is written to a .mat file in `labels_data`.

It also moves the matching .lst file to `training_test_data`.

### Files, Directories

Directory `training_test_data` contains pedobarograph (pressure pad) data.

- Each `.lst` file represents one reading of a single footstep on the pedobarograph. 
It is read using `pedo_extract.m`.

Directory `labels_data` contains labeled matrices.

- Each `.mat` file is a labelled matrix corresponding to the max-pressure map 
of the `.lst` file with the same id, eg `pb_data/NP10_2651.lst` corresponds to `labels_data/NP10_2651.mat`.
- The labelled matrix contains 1s for the great toe, 2s for the lateral forefoot, 
3s for the medial forefoot, 4s for the heel, and 0s for all other entries.

File `pedo_extract.m` (Written by: Quinn Boser, July 2013) takes a `lst` file and returns

1. a 3d struct array with dimensions (x,y,time)
2. a 2d max-pressure map with dimensions (x,y)
3. the length of the x dimension
4. the length of the y dimension
5. the length of the time dimension

File `pedo_format.txt` describes the format of an `.lst` file. It is used by `pedo_extract.m`.

File `get_params.m` reads the number of rows, columns, and frames from an `.lst` file, where frames is the length of the time dimension. It is used by `pedo_extract.m`.

File `JaccardScore.m` takes 2 matrices, A and B, and returns their Jaccard Score.

File `silly_example.m` is a small example of training a Conditional Random Field (CRF) on noise. It is pulled from http://users.cecs.anu.edu.au/~jdomke/JGMT/#binarydenoising under the MIT License.

File `sophisicated_example.m` is an example of training a Conditional Random Field (CRF) on the Stanford Backgrounds dataset. It is pulled from http://users.cecs.anu.edu.au/~jdomke/JGMT/#backgrounds under the MIT License.

