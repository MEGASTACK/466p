function [pedo_dat,pedo_max,rows,columns,frames] = pedo_extract(filename )
% pedo_extract: Determines parameters needed to extract information from
%              pedobarograph data file (number of frames, size of
%              matrices). Extracts data from all time frames in the data
%              file and generates a single matrix representing the maximum
%              pressure registered on each cell throughout the step.
%              Displays maximum pressure matrix (foot image).
% Written by: Quinn Boser, July 2013
%   
%   Syntax: [pedo_dat,pedo_max,rows,columns,frames] = pedo_extract( filename, handles )
%
%   Input:  filename = string with name of .lst pedobarograph ascii data
%                      file
%           ***REMOVED *** handles = structure with handles and user data from the pedo
%           GUI
%
%   Output: pedo_dat = cell array of data sets from each time frame in file 
%           pedo_max = maximum pressure matrix
%           rows = integer number of rows in each time frame matrix and in
%                  pedo_max
%           columns = integer number of columns in each time frame matrix
%                     and in pedo_max
%           frames = number of time frames collected in data file
%   
%   Functions and files required to run: 
%           get_params.m
%           pedo_format.txt
%           custom color map pedcmap1
%
%   Required toolboxes:
%           Image Processing Toolbox

% cla(handles.maxPfoot)                               % clear handles on GUI

%---------------------------------------------------% 
% Open data file and determine paremeters:
%---------------------------------------------------% 
pedoID = fopen(filename);                           % open pedo-barograph data file  
[rows,columns,frames,up_down] = get_params(pedoID);

%---------------------------------------------------% 
% Read row format based on # of columns:
%---------------------------------------------------%
formatID = fopen('pedo_format.txt');                % load pedo_format.txt
for k=1:columns                                     % use number of columns in matrix to...
    formatscan = fgetl(formatID);                   % ...obtain format of rows
end
fclose(formatID);                                   % close pedo_format.txt

%---------------------------------------------------% 
% Allocate memory:
%---------------------------------------------------% 
pedo_dat = cell(1,frames);                          % cell array for all data                              
pedo_max = zeros(rows,columns);                     % matrix of maximum pressures

%---------------------------------------------------% 
% Extract data from file:
%---------------------------------------------------% 
 for n = 1:frames                                   % run through every frame (time instance) in data file
     pedo_cell = textscan(pedoID,formatscan,rows,'HeaderLines',12);     % store current frame from data file (will be in format cell array)
      if (up_down)                                  % if image will be up-side-down...
        pedo_new = rot90(cell2mat(pedo_cell),2);    % ...convert cell array to matrix and rotate 180 degrees
      else                                          % if image will right side up...
        pedo_new = cell2mat(pedo_cell);             % ...convert cell array to matrix 
      end
     for i = 1:rows                                 % run through every cell in matrix
         for j = 1:columns
             if (pedo_new(i,j) > pedo_max(i,j))     % if value in cell is greater than previous maximum
                 pedo_max(i,j) = pedo_new(i,j);     % store it as new maximum
             end
         end
     end        
     pedo_dat{n} = pedo_new;                        % store matrix in cell array (Note: not entirely necessary if memory becomes issue)
     textscan(pedoID,'%*[^\n]',2);                  % skip last line of data
 end
 
%---------------------------------------------------% 
% Generate maximum pressure image:
%---------------------------------------------------% 
% axes(handles.maxPfoot);
% imagesc(pedo_max,[0,300])                                   % plot maximum matrix 
% set(handles.maxPfoot,'XLim',[0.5,columns+0.5],'YLim',[0.5,rows+0.5])
% hold on
% load('MyColormaps','pedcmap3');                     % load custom color map: pedomap1
% colormap(pedcmap3);                                 % apply color map to figure
% colorbar;                                           % view color bar
% for i = 1:rows                                      % run through all cells in maximum pressure matrix
%     for j = 1:columns
%         if pedo_max(i,j) ~= 0                       % if pressure value is not zero
%             pedo_str = num2str(pedo_max(i,j));      % turn value into string
%             text(j,i,pedo_str,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',7);    % write value over cell on figure
%         end
%     end
% end

fclose(pedoID);                                     % close data file

end

