# T10-Gait

Foot segment identification from footsteps on pressure pads

Google Drive folder is here: https://drive.google.com/open?id=0B2WhXU_IjE68cjlmMzZZdV9iMm8

Proposal is here: https://drive.google.com/open?id=1vuIcVz8KHbXqyjY4PWqkA_2YfWaEg1bcKuxf7oa9vcg

### Usage

File `pedo_extract.m` takes a `lst` file and returns:
1. a 3d struct array with dimensions (x,y,time)
2. a 2d max-pressure map with dimensions (x,y)

