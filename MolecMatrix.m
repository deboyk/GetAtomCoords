function [ M, Atoms, numatoms, Connectivity] = MolecMatrix( filename )
%MolecMatrix Summary: Function to calculate the molecular matrix (M)
% This function extracts coordinates from a .sdf file and generates the
% molecular matrix (M). This function also extracts the connectivity table
% and creates a separate matrix for this information. The connectivity
% table contains the bond information.
%
% See Also: R3mCalculate, ImportSDF
%
% Author: Kevin DeBoyace
%         Wildfong Lab
%         Duquesne University
% Updated: Jan 2019

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
delimiter = ' ';
startRow = 2;
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'ReturnOnError', false);
fclose(fileID);
for ii = 1:size(dataArray{1,1},1)
    for jj = 1:size(dataArray,2)
        dataArraynum(ii,jj) = str2double(dataArray{1,jj}{ii});
    end
end

%% Extract cartesian coordinates from file and place them in a single matrix.
dataArraynum_cut = dataArraynum(:,1:3); %Take first three columns
datarem = rem(dataArraynum_cut,1);
jj = 1;
for ii = 1:size(dataArraynum_cut, 1);
    if sum(datarem(ii,:)) == 0 || isnan(datarem(ii,1))
        temp_rem(ii) = 0; 
    else
        temp_rem(ii) = 1; 
        coords(jj,:) = dataArraynum_cut(ii,:);
        Atoms(jj,:) = dataArray{1,4}(ii); % Extract atom names
    end
    jj = jj+1;
end

for ii = 1:size(Atoms,1)
    temp_emp(ii) = isempty(Atoms{ii});
    try
        temp(ii) = str2num(Atoms{ii});
    catch
        temp(ii) = 0;
    end
end
[row,col] = find(temp ~= 0 | temp_emp == 1);
Atoms(col) = [];
coords(col,:) = [];

M = coords; % Matrix of cartesian coordinates
clear ii jj

numatoms = size(M,1);

%% Connectivity table
jj = 1;

for ii = 1:size(dataArraynum_cut,1)
    if (sum(datarem(ii,:))) == 0 && ~isnan(datarem(ii,1)) && sum(dataArraynum_cut(ii,1:3)~=0) == 3 %For those numbers which are part of connectivity or NaN --> delete
        e(jj,:) = dataArraynum_cut(ii,:); %Build matrix with connectivity data
        jj = jj + 1;
    else
        jj = jj;
    end
end

Connectivity = e;

clear ii jj d e f %Clear Temp variables

end

