function [ mp_pixel ] = midpoint_pixel( input_arr, i, max_val )
%MIDPOINT_PIXEL Returns the integer midpoint between input_arr(i) and input_arr(i+1)
% Will not return a value larger than max_val
    mp_pixel = floor(input_arr(i) + (input_arr(i+1) - input_arr(i)));
    if mp_pixel > max_val
        mp_pixel = max_val;
    end
end

