function label_foot( np_and_number, foot_number )%, output_filename )
%LABEL_FOOT
    
    addpath(['..' filesep]); % for pedo_extract
    
    output_path = ['..' filesep 'labels_data' filesep];
    input_path = ['..' filesep 'FinalData' filesep np_and_number filesep foot_number];
    imgpath = [input_path '.jpg'];
    lstpath = [input_path '.lst'];
    
    img = imread(imgpath);
    [lst, lst_max, lst_rows, lst_cols, lst_frames] = pedo_extract(lstpath);
    
    disp_step(''); % reset step counter
    
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
    
    regions_mask = zeros(size(img,1), size(img,2));
    
    disp_step('Trace the toe region (region 1).');
    disp('Draw a closed loop to continue.');
    toe_mask = roipoly(img);
    regions_mask(toe_mask==1) = 1;
    close all;
    
    disp_step('Trace the medial forefoot region (region 2). It is below the toe.');
    disp('Draw a closed loop to continue.');
    m_forefoot_mask = roipoly(img);
    regions_mask(m_forefoot_mask==1) = 2;
    close all;
    
    disp_step('Trace the lateral forefoot region (region 3). It is beside the medial forefoot.');
    disp('Draw a closed loop to continue.');
    l_forefoot_mask = roipoly(img);
    regions_mask(l_forefoot_mask==1) = 3;
    close all;
    
    disp_step('Trace the heel region (region 4). It is at the bottom.');
    disp('Draw a closed loop to continue.');
    heel_mask = roipoly(img);
    regions_mask(heel_mask==1) = 4;
    close all;
    
    foot_right = foot_left + foot_width;
    foot_bottom = foot_top + foot_height;
    block_width = foot_width / lst_cols;
    half_block_width = block_width / 2;
    block_height = foot_height / lst_rows;
    half_block_height = block_height / 2;
    
    out_mat = zeros(lst_rows, lst_cols);
    
    foot_pixel_collist = foot_left:block_width:foot_right;
    foot_pixel_rowlist = foot_top:block_height:foot_bottom;
    lst_rowlist = 1:lst_rows;
    lst_collist = 1:lst_cols;
    for i=lst_rowlist
        mid_point_row = midpoint_pixel(foot_pixel_rowlist, i, size(regions_mask, 1));
        for j=lst_collist
            mid_point_col = midpoint_pixel(foot_pixel_collist, j, size(regions_mask, 2));
            out_mat(i,j) = regions_mask(mid_point_row, mid_point_col);
        end
    end
    
    save([output_path np_and_number '_' foot_number '.mat'], 'out_mat');
end
