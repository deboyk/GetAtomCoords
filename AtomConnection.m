function [ connections, topology_mat ] = AtomConnection( Connectivity, Atoms )
%AtomConnection - A function which determines which atoms are connected by
%distances of 1, 2 and 3 bond distances.
%   This function takes in the connectivity table (without the bond order
%   column), and outputs a structure containing cell arrays, where each
%   column corresponds to the atom number, and each row corresponds to the
%   number of bond distances (e.g. col 2, row 3 represents the atoms which
%   are 3 bond distance away from atom 2.) 
%
%   See also: ImportSDF, PlotMolecule, R3mCalculate
%
%   Author: Kevin DeBoyace & Shikhar Mohan
%   Updated: Jan 2019


connect = Connectivity(:,1:2);

for i = 1:size(Atoms,1)
    B = find(connect(:,1) == i);
    C = find(connect(:,2) == i);
    Temp.A{1,i} = [connect(B,:);connect(C,:)];
    Temp.A{1,i} = reshape(Temp.A{1,i},1,[]);
    for j = 1:size(Temp.A{1,i},2);
        if Temp.A{1,i}(j) == i
            Temp.A{1,i}(j) = 0; %Replace where row = value with zeros (Itself present in its own row)
        else
        end
    end
    %Delete zeros 
    test = Temp.A{1,i};
    test(test==0) = [];
    Temp.A{1,i} = test;
end
% remove empty cells

clear B C i j test


% Working 2 bond connection - Need to delete atom self reference

for i=1:size(Atoms,1)
    for jj = 1:size(Temp.A{1,i},2)
        B = find(connect(:,1) == Temp.A{1,i}(jj));
        C = find(connect(:,2) == Temp.A{1,i}(jj));
         for kk = 1:size(connect(B,:),2) 
                D{i, jj} = reshape([connect(B,:);connect(C,:)], 1, []);
        end
    end
end

% Place values into Temp structure
for i = 1:size(D,1);
    Temp.A{2,i} = [D{i,:}];
end
clear B C D E F i ii j jj k kk 

% Delete atoms which are 1 bond distance away.
for i = 1:size(Temp.A,2)
    for k = 1:size(Temp.A{1,i},2)
        [~, col] = find(Temp.A{1,i}(k) == Temp.A{2,i});
        Temp.A{2,i}(col) = [];
    end
end
clear i k col row

% Delete where column = atom number (where atom itself is included in bond
% connections)

for i = 1:size(Temp.A,2)
    if isempty(Temp.A{2,i}) == 0; % NEW
        a = find(Temp.A{2,i} == i);
        Temp.A{2,i}(a) = [];
    elseif isempty(Temp.A{2,i}) == 1; % NEW
	    Temp.A{2,i}(a) = []; % NEW

end


%  3 Bond Distances 

for i=1:size(Atoms,1) % NEW
    for jj = 1:size(Temp.A{2,i},2)
        if isempty(Temp.A{2,i}) == 0;
            B = find(connect(:,1) == Temp.A{2,i}(jj));
            C = find(connect(:,2) == Temp.A{2,i}(jj));
             for kk = 1:size(connect(B,:),2) 
                if connect(B,kk) ~= i & connect(B,kk) ~= jj 
                    D{i, jj} = reshape([connect(B,:);connect(C,:)], 1, []);
                else
                end                
             end
        else
        end        
    end
end


% Place values into Temp structure
for i = 1:size(D,1);
    Temp.A{3,i} = [D{i,:}];
end
clear B C D E F i ii j jj k kk a


% Delete atoms which are 1 bond distance away.
for i = 1:size(Temp.A,2)
    for k = 1:size(Temp.A{1,i},2)
        [~, col] = find(Temp.A{1,i}(k) == Temp.A{3,i});
        Temp.A{3,i}(col) = [];
    end
end
clear i k col row


% Delete atoms which are 2 bond distance away.
for i = 1:size(Temp.A,2)
    for k = 1:size(Temp.A{2,i},2)
        [~, col] = find(Temp.A{2,i}(k) == Temp.A{3,i});
        Temp.A{3,i}(col) = [];
        
    end
end
clear i k col row

% Delete duplicate values and sort values
for ii = 1:size(Temp.A,1)
    for jj = 1:size(Temp.A, 2)
        Temp.A{ii,jj} = unique(Temp.A{ii,jj}); % Keeps only unique values.
    end
end

connections = Temp;


% Topological distance matrix
topology_mat = zeros(size(Connectivity,1)); %Preallocate matrix

    for ii = 1:size(Atoms,1)
        for jj = 1:size(connections.A{3,ii},2)
            %try
            if isempty(connections.A{3,ii}) == 0; 
                topology_mat(ii,connections.A{3,ii}(jj)) = 1;
                topology_mat(connections.A{3,ii}(jj), ii) = 1;
            else 
            end
        end
    end



end




