# T10-Gait

Foot segment identification from footsteps on pressure pads

Google Drive folder is here: https://drive.google.com/open?id=0B2WhXU_IjE68cjlmMzZZdV9iMm8

Proposal is here: https://drive.google.com/open?id=1vuIcVz8KHbXqyjY4PWqkA_2YfWaEg1bcKuxf7oa9vcg

### Usage

File `main.m`

### Files, Directories

Directory `pb_data` contains pedobarograph (pressure pad) data.
- Each `lst` file represents one reading of a single footstep on the pedobarograph. It is read using `pedo_extract.m`.

File `pedo_extract.m` (Written by: Quinn Boser, July 2013) takes a `lst` file and returns
1. a 3d struct array with dimensions (x,y,time)
2. a 2d max-pressure map with dimensions (x,y)
3. the length of the x dimension
4. the length of the y dimension
5. the length of the time dimension

File `JaccardScore.m` takes 2 matrices, A and B, and returns their Jaccard Score.

File `silly_example.m` is a small example of training a Conditional Random Field (CRF) on noise. It is pulled from http://users.cecs.anu.edu.au/~jdomke/JGMT/#binarydenoising under the MIT License.

File `sophisicated_example.m` is an example of training a Conditional Random Field (CRF) on the Stanford Backgrounds dataset. It is pulled from http://users.cecs.anu.edu.au/~jdomke/JGMT/#backgrounds under the MIT License.

