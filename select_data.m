function [ images_names, labels_names ] = select_data( images_path, labels_path, options)
%SELECT_DATA returns lists of data to be processed.
%   Detailed explanation goes here

all = 'all';
left = 'left';
right = 'right';
one = 'one';


all_images = dir([images_path '*.lst']);
all_labels = dir([labels_path '*.mat']);
images_names = [];
labels_names = [];


if (~exist('options', 'var'))
    options = all;
end



if (strcmp(options, all))
    images_names = all_images;
    labels_names = all_labels;
elseif (strcmp(options, one))
    
    found_im = [];
    found_lab = [];
    
    for im = all_images'
        filename = im.name;
        if(length(strfind(found_im, filename(1:4))) == 0)
            found_im = [found_im filename(1:4)];
            images_names = [images_names, im];
        end
    end
    
    for lab = all_labels'
        filename = lab.name;
        if(length(strfind(found_lab, filename(1:4))) == 0)
            found_lab = [found_lab filename(1:4)];
            labels_names = [labels_names, lab];
        end
    end
    
    
    
end



end

