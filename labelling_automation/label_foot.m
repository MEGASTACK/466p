function label_foot( np_and_number, foot_number )%, output_filename )
%LABEL_FOOT
    
    %% 
    % Includes
    %
    addpath(['..' filesep]); % for pedo_extract
    
    %%
    % Constants
    %
    PROJ_ROOT = ['..' filesep];
    LABELS_DIR = [PROJ_ROOT 'labels_data' filesep];
    FINAL_DATA_DIR = [PROJ_ROOT 'FinalData' filesep];
    TRAIN_TEST_DIR = [PROJ_ROOT 'train_test_data' filesep];
    
    
    %%
    % Read input, build paths
    %
    output_path = LABELS_DIR;
    input_path =  [FINAL_DATA_DIR np_and_number filesep foot_number];
    imgpath = [input_path '.jpg'];
    lstpath = [input_path '.lst'];
    
    %%
    % Error checking
    %
    error_should_exit = false;
    if exist(imgpath, 'file')
        fprintf(['Reading labelled image ' imgpath '.\n']);
    else
        fprintf(['Labelled image ' imgpath ' does not exist\n']);
        error_should_exit = true;
    end
    if exist(lstpath, 'file')
        fprintf(['Reading footstep file ' lstpath '.\n']);
    else
        fprintf(['Footstep file ' lstpath ' does not exist\n']);
        error_should_exit = true;
    end
    if error_should_exit
        return
    end
    fprintf(['\n']); % separate IO output from rest of program outpu
    
    %%
    % Main Functionality begins
    %
    
    %%
    % Read the jpg (image), read the lst (footstep)
    %
    img = imread(imgpath);
    [lst, lst_max, lst_rows, lst_cols, lst_frames] = pedo_extract(lstpath);
    
    disp_step(''); % reset step counter
    
    %%
    % Get the bounding box of the foot. The lst file is bounded to the
    % region which contains all non-zero squares in the max-pressure matrix.
    %
    % So, this step instructs the user to move and resize a rectangle to match 
    % the region which contains all non-zero squares in the max-pressure matrix.
    %
    disp_step('Adjust the rectangle to match the outline of the colored squares.');
    disp('Double-click inside the rectangle to continue.');
    figure, imshow(img, 'InitialMagnification', 150);
    h = imrect(gca, [40 40 size(img,2)-40 size(img,1)-40]);
    foot_box = wait(h);
    foot_left = foot_box(1);
    foot_top = foot_box(2);
    foot_width = foot_box(3);
    foot_height = foot_box(4);
    close all;
    
    %%
    % This regions mask will be mapped later to the dimensions of the lst
    % file.
    %
    % It contains the following regions, marked by the following values:
    % - [0] N/A
    % - [1] Toe
    % - [2] Medial Forefoot
    % - [3] Lateral Forefoot
    % 
    % Each entry in this mask corresponds to a pixel in the jpg image.
    %
    regions_mask = zeros(size(img,1), size(img,2));
    
    %%
    % Prompt the user to trace the toe region. This uses MATLAB builtin
    % function `roipoly` which opens the image, allows the user to click to
    % draw a polygon, then returns a logical matrix of the same dimensions
    % as the image. The matrix contains 1 for (x,y) coordinates which are
    % *inside* the drawn region, and 0 for (x,y) coordinates which are
    % *outside* the drawn region.
    %
    disp_step('Trace the toe region (region 1).');
    disp('Draw a closed loop to continue.');
    toe_mask = roipoly(img);
    regions_mask(toe_mask==1) = 1;
    close all;
    
    %%
    % Prompt the user to trace the medial forefoot.
    % Also describe where the medial forefoot is.
    %
    disp_step('Trace the medial forefoot region (region 2). It is below the toe.');
    disp('Draw a closed loop to continue.');
    m_forefoot_mask = roipoly(img);
    regions_mask(m_forefoot_mask==1) = 2;
    close all;
    
    %%
    % Prompt the user to trace the lateral forefoot.
    % Also describe where the lateral forefoot is.
    %
    disp_step('Trace the lateral forefoot region (region 3). It is beside the medial forefoot.');
    disp('Draw a closed loop to continue.');
    l_forefoot_mask = roipoly(img);
    regions_mask(l_forefoot_mask==1) = 3;
    close all;
    
    %%
    % Prompt the user to trace the heel.
    %
    disp_step('Trace the heel region (region 4). It is at the bottom.');
    disp('Draw a closed loop to continue.');
    heel_mask = roipoly(img);
    regions_mask(heel_mask==1) = 4;
    close all;
    
    % Calculate some useful values in preparation for mapping the
    % regions_mask into a matrix with the dimensions of the lst file
    foot_right = foot_left + foot_width;
    foot_bottom = foot_top + foot_height;
    block_width = foot_width / lst_cols;
    half_block_width = block_width / 2;
    block_height = foot_height / lst_rows;
    half_block_height = block_height / 2;
    
    % initialize the output file
    out_mat = zeros(lst_rows, lst_cols);
    
    %%
    % foot_pixel_*list arrays are 1 longer than their lst_*list
    % counterparts. This is because each lst entry corresponds to the
    % *center* of two foot_pixel entries.
    %
    % eg
    %
    % |o|o|o|o|
    % 
    % o lst
    % | foot_pixel
    %           
    foot_pixel_collist = foot_left:block_width:foot_right;
    foot_pixel_rowlist = foot_top:block_height:foot_bottom;
    lst_rowlist = 1:lst_rows;
    lst_collist = 1:lst_cols;
    
    %%
    % Error checking: assert that foot_pixel_lists and lst_lists are of the
    % correct size. If they aren't either something is wrong with the lst
    % file or the user's bounding rectangle was not good.
    %
    if (length(lst_collist)+1 ~= length(foot_pixel_collist))
        fprintf('Footstep file dimensions and bounding rectangle dimensions don''t match\n');
        fprintf('Exiting.\n');
    end
    
    %%
    % For each midpoint between foot_pixel values, get the pixel in the
    % middle and save it to the corresponding index in the output matrix.
    % This maps the regions_mask from pixel-dimensions into lst-dimensions.
    %
    for i=lst_rowlist
        mid_point_row = midpoint_pixel(foot_pixel_rowlist, i, size(regions_mask, 1));
        for j=lst_collist
            mid_point_col = midpoint_pixel(foot_pixel_collist, j, size(regions_mask, 2));
            out_mat(i,j) = regions_mask(mid_point_row, mid_point_col);
        end
    end
    
    %%
    % Write the labelled foot data to file
    %
    % Filenames are of format /output/path/NP##_####.mat
    %
    fprintf('\n'); % separate IO output from rest of program
    outfilename = [output_path np_and_number '_' foot_number '.mat'];
    save(outfilename, 'out_mat');
    fprintf(['Success! Saved labelled foot as ' outfilename '\n']);
    
    %%
    % Move the consumed lst file to the training_test_data directory. Now that it
    % is labelled, it is available to be trained with or tested on.
    %
    lstmovepath = [TRAIN_TEST_DIR np_and_number '_' foot_number '.lst'];
    movefile(lstpath, lstmovepath);
    fprintf(['Also moved footstep from ' lstpath ' to ' lstmovepath '\n']);
end
