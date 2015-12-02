function [rows,columns,frames,up_down] = get_params(fileID)
% get_params : determines parameters needed to extract data from
% pedo-barograph ascii file
% Written by: Quinn Boser, June 2013
%
%   Syntax: [rows,columns,frames,up_down] = get_params(fileID)
%
%   Input:  fileID = integer that identifies the pedobarograph data file
%                   (output of fopen)
%         
%   Output: rows = number of rows in each time instance matrix in the pedo
%                   barograph data file
%           columns = number of columns in each matrix in file
%           frames = number of data sets (time instances) collected
%           up_down = 0 if matrices are oriented with heel data at the
%                     bottom
%                   = 1 if matrices are oriented with heel data at the top 

for k = 1:12
    column_label = fgetl(fileID);                   % obtain 12th line of file
end
columns = (length(column_label)-13)/8;              % determine # of columns in foot matrix

check = 1;                                          % initiate check as 1
rows = 0;                                           % initiate # of rows at 0
found_start_side = 0;
while (check==1)                                    % read lines until first token in line = "Force":
    row = fgetl(fileID);                            % read row 
    [first,rest_row] = strtok(row);                 % split first token from rest
    if (strcmp(first,'Force'))                      % if first token is "Force"
        check = 0;                                  % stop reading rows
    else
        rows = rows+1;                              % increment # of rows if haven't found "Force" yet
        rest_arr = str2num(rest_row);               % convert rest of row to number
        if (max(rest_arr) && found_start_side==0)   % if first row with non-zero value
            start_side = rows;                      % store row number
            found_start_side = 1;                   % indicate that first row with data has been found
        end     
    end
end

if (start_side > rows/2)                            % if heel at bottom of matrix...
    up_down = 0;                                    % ...image will be right-side-up (set up_down = 0)
else                                                % if heel at top of matrix...
    up_down = 1;                                    % ...image will be up-side-down (set up_down = 1)
end

rest_text = textscan(fileID,'%s','delimiter','\n'); % read rest of file
lines = length(rest_text{1});                       % determine # lines in file
frames = floor(lines/(rows+13))+1;                  % determine # time frames in file

frewind(fileID);                                    % set pointer to beginning of data file


end

