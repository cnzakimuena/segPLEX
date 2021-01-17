
function call_segPLEX()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name - call_segPLEX()
% Creation Date - 15th January 2021
% Author - Charles Belanger Nzakimuena
% Website - https://www.ibis-space.com/
%
% Description - 
%   PRESEG performs initial segmentation prior to refined segmentation
%   folderPath : path towards patient folders.
%   All the data is saved in each patient folder in the subfolder 'DataFiles'
%
%   REFSEG performs refined segmentation prior to network graph implementation
%   folderPath : path towards patient folders.
%   All the data is saved in each patient folder in the subfolder 'Results'
%
%   ETDRSGRIDS generates 2D and 3D ETDRS grids
%
% addpath(folderName1,...,folderNameN) adds the specified folders to the 
% top of the search path for the current MATLAB session
% genpath(folderName) returns a character vector containing a path name 
% that includes folderName and multiple levels of subfolders below 
% folderName (relative directory indicated by the dot symbol is collapsed 
% with the the specified folder, i.e., the subfolder 'subfunction')
%
% Example -
%		call_segPLEX()
%
% License - MIT
%
% Change History -
%                   15th January 2020 - Creation by Charles Belanger Nzakimuena
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('./subfunctions'))

%% list names of folders inside the patients folder
 
currentFolder = pwd;
patients_Folder = fullfile(currentFolder, 'preprocessed');
preSeg(patients_Folder)
refSeg(patients_Folder)
ETDRSgrids(patients_Folder)



