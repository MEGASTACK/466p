function [ filenames ] = GetLstFilenames( path )
%% GetLstFiles
%   Given the root directory <path>, return a cell containing the path to
%   all .lst files in that directory and all sub-directories.

root = dir(path);
numFiles = 1;
fNames = cell(0, 0);
recursiveFNames = cell(0, 0);

for i = 1 : length(root)
    if or(strcmp(root(i).name, '.'), strcmp(root(i).name, '..'))
        continue;
    end
    
    [fDir, fName, ext] = fileparts([path, '/', root(i).name]);
    
    if and(root(i).isdir, fName)
        recursiveFNames = [recursiveFNames; GetLstFilenames([fDir, '/', fName])];
    end
    
    if strcmp(ext, '.lst')
        fNames{numFiles, 1} = [path, '/', fName, ext];
        numFiles = numFiles + 1;
    end
end

filenames = [fNames; recursiveFNames];

end
