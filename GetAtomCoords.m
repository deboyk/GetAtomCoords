function [r3m_out] = GetAtomCoords()
%GetAtomCoords_multi - This function allows the user to select data exported from
%multiple Materials Studio files (.xsd files) and calculate R3m values
%for each molecule in each frame. A distribution of R3m values is then
%plotted as histogram.
%   This function allows the user to select text files to
%   determine the distribution of R3m values that result from
%   conformational differences in molecular structures as determined from
%   Molecular Dynamics Simulations performed using Materials Studio. 
%
%   See also: R3mCalculate_auto, AtomConnection, AtomicWeighting,
%   EuclidDistance,InfluenceDistanceMat, MolecInfluenceMatrix, MolecMatrix
%
%   Code Author: Kevin DeBoyace
%                Duquesne University
%   
%   Last Updated:   January 2019
%   
%   NOTE: The Materials Studio data must be saved as individual 3D 
%   atomistic files (.xsd format).
% 
% 
%% IMPORT MULTIPLE FILES
[files_in, path_in] = uigetfile('*.xsd','Select the INPUT DATA FILE(s)', 'Multiselect', 'on');

%% Wait bar
w = waitbar(0,'Initiating...');

%% FOR 1 FILE SELECTED
if isstr(files_in) == 1;

    file_char = ([path_in files_in]); % Import multiple files
    filename = file_char;
    filecontent = fileread(filename);

    %% Wait bar
    w = waitbar(0,'Wait');

    %% Set number of molecules in system
    num_mol = 40;

    %% Get Data
    data1 = regexp(filecontent, '<Atom3d ID=...........................................................................................................................................................................................................................................', 'Match')';
    data2 = regexp(filecontent, '<Bond ID=.................................................................', 'Match')';

    % Atom ID
    atom_id = extractBetween(data1, '<Atom3d ID="', '"');
    atom_id = str2double(atom_id); % convert to numeric

    % Molecule ID
    molecule_id = extractBetween(data1, ' Parent="', '"');
    molecule_id = str2double(molecule_id); % convert to numeric
    
    min_mol = min(molecule_id);
    max_mol = max(molecule_id);

    % Coordinates
    xyz = extractBetween(data1, ' XYZ="', '"');
    coords = regexp(xyz, ',', 'split');
    coords = vertcat(coords{:});
    coords = str2double(coords); % convert to numeric
    clear xyz

    % Get Box limits
    max_x = max(coords(:,1));
    min_x = min(coords(:,1));
    max_y = max(coords(:,2));
    min_y = min(coords(:,2));
    max_z = max(coords(:,3));
    min_z = min(coords(:,3));

    % Get Cell Parameters
    data3 = regexp(filecontent, 'AVector=...............................', 'Match')';
    cell_param = extractBetween(data3, '"', ',0,0"');
    cell_param = str2double(cell_param);

    % Get Bond ID
    forBond = regexp(filecontent, '<Bond ID=..........', 'Match')';
    bond_id = extractBetween(forBond, '<Bond ID="', '"');
    clear forBond

    % Get Connections
    connections = extractBetween(data2, 'Connects="','"');
    connect = regexp(connections, ',', 'split'); %column 1 = x, column 2 = y, column 3 = z
    connect = vertcat(connect{:});
    connect = str2double(connect); % convert to numeric
    num_connect = size(connect,1)/num_mol;

    % Get Atoms
    forAtom = regexp(filecontent, 'Components=........', 'Match')';
    atom = extractBetween(forAtom, 'Components="','"');
    clear forAtom

    connect_actual = 1:num_connect;

    minMol = min(molecule_id);
    maxMol = max(molecule_id);
    allMol = unique(molecule_id); % identify all unique values in 'molecule_id'

    % Save to .xyz format
    %for hh = minMol:maxMol;
    for hh = 1:size(allMol,1) 
        mol1 = find(molecule_id == allMol(hh));
        waitbar(hh/maxMol, w, 'Calculating...');
        %mol1 = find(molecule_id == hh);
        atom_id_1 = atom_id(mol1);
        coords_1 = coords(mol1,:);
        connect_1 = connect(connect_actual,:);
        connect_actual = connect_actual(end):connect_actual(end)+num_connect;
        atom_1 = atom(mol1);
        for ii = 1:size(atom_id_1,1);
            % ii = atom row of interest
            [row, ~] = find(atom_id_1(ii) == connect_1); %  row in connect where atom is found (column 1)
            connected = connect_1(row,:);
            % Get coordintes for relevant atoms
            for jj = 1:size(connected,1);
                for kk = 1:size(connected,2);
                    coord_find(jj,kk) = {coords_1(find(atom_id_1==connected(jj,kk)),:)};
                end
            end
            % XYZ format
            % number of elements
            % comment line
            % <elemente> <X> <Y> <Z>
            % ...

            % First line: num_atoms;
            % comment line
            % atom_1; coords_1;
            num_atoms = size(atom_1,1);
            atom_label = char(atom_1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % NOTE: data normalized to 1 (i.e. relative distance is the same, but overall
            % distance is no longer representative. Must convert back to Angstroms
            % HOW? ---> USE CELL VOLUME
            coords_conv = coords_1*cell_param; % Use cell parameter to convert back to angstroms
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            file_xyz = ['testOut',num2str(hh),'.xyz'];
            fileID = fopen(file_xyz, 'wt');
            fprintf(fileID, '%d\n', num_atoms);
            fprintf(fileID, 'NAME\n');
            for ii = 1:size(atom_1,1);
                %fprintf(fileID, '%s\t', '%i\t %i\t% i\t\n', atom_label, coords_1);
                fprintf(fileID, '%s\t %i\t %i\t% i\t\n',atom_label(ii,:), coords_conv(ii,:));
            end
            fclose(fileID);
            clear fileID
        end
        %% Convert from xyz to sdf to allow calculation of R3m
        % Run Open Babel ot convert
        % see https://openbabel.org/docs/dev/Command-line_tools/babel.html
        % -c = center coordinates
        % -o = specify output format (e.g. -osdf --> output to sdf)
        copyfile(['testOut',num2str(hh),'.xyz'], 'testOut.xyz'); %Copy file to run in Babel
        %!obabel testOut.xyz -c -O Out.sdf
        [yy, zz] = system('obabel testOut.xyz -c -O Out.sdf');
        clear yy zz
        copyfile('Out.sdf', ['testOut',num2str(hh),'.sdf']); % Copy output file from Babel to new file name
        %% Calculate R3m
        pathname = cd;
        file_sdf = ['testOut',num2str(hh),'.sdf'];
        filename = [pathname,'\', file_sdf]; % stitch together for full path
        [r3m] = R3mCalculate_auto(filename);
        %r3m_out((hh-minMol+1),1) = r3m; % Output r3m value to vector
        r3m_out(hh,1) = r3m;
    
end

%% FOR MULTIPLE FILES SELECTED 
else
    files_in = files_in';
    for aa = 1:size(files_in,1);

        file_char = ([path_in files_in{aa}]); % Import multiple files
        filename = file_char;
        filecontent = fileread(filename);

        % Set number of molecules in system
        num_mol = 40;

        % Get Data
        data1 = regexp(filecontent, '<Atom3d ID=...........................................................................................................................................................................................................................................', 'Match')';
        data2 = regexp(filecontent, '<Bond ID=.................................................................', 'Match')';

        % Get Atom ID
        atom_id = extractBetween(data1, '<Atom3d ID="', '"');
        atom_id = str2double(atom_id); % convert to numeric

        % Get Molecule ID
        molecule_id = extractBetween(data1, ' Parent="', '"');
        molecule_id = str2double(molecule_id); % convert to numeric

        % Get Coordinates
        xyz = extractBetween(data1, ' XYZ="', '"');
        coords = regexp(xyz, ',', 'split');
        coords = vertcat(coords{:});
        coords = str2double(coords); % convert to numeric
        clear xyz

        % Get Box limits
        max_x = max(coords(:,1));
        min_x = min(coords(:,1));
        max_y = max(coords(:,2));
        min_y = min(coords(:,2));
        max_z = max(coords(:,3));
        min_z = min(coords(:,3));

        % Get Cell Parameters
        data3 = regexp(filecontent, 'AVector=...............................', 'Match')';
        cell_param = extractBetween(data3, '"', ',0,0"');
        cell_param = str2double(cell_param);

        % Get Bond ID
        forBond = regexp(filecontent, '<Bond ID=..........', 'Match')';
        bond_id = extractBetween(forBond, '<Bond ID="', '"');
        clear forBond

        % Get Connections
        connections = extractBetween(data2, 'Connects="','"');
        connect = regexp(connections, ',', 'split'); %column 1 = x, column 2 = y, column 3 = z
        connect = vertcat(connect{:});
        connect = str2double(connect); % convert to numeric
        num_connect = size(connect,1)/num_mol;

        %NOTE: Connections appear to refer to bond_id    

        % Get Atoms
        forAtom = regexp(filecontent, 'Components=........', 'Match')';
        atom = extractBetween(forAtom, 'Components="','"');
        clear forAtom

        % Select only relevent connections (b/c numatoms ~= num connections)
        connect_actual = 1:num_connect;

        % Identify min and max molecule id numbers for iteration
        minMol = min(molecule_id);
        maxMol = max(molecule_id);
        allMol = unique(molecule_id); % identify all unique values in 'molecule_id'

        %% Save to .xyz format
        for hh = 1:size(allMol,1) 
            mol1 = find(molecule_id == allMol(hh));
            atom_id_1 = atom_id(mol1);
            coords_1 = coords(mol1,:);
            connect_1 = connect(connect_actual,:);
            connect_actual = connect_actual(end):connect_actual(end)+num_connect;
            atom_1 = atom(mol1);

            for ii = 1:size(atom_id_1,1);
                [row, ~] = find(atom_id_1(ii) == connect_1); %  row in connect where atom is found (column 1)
                connected = connect_1(row,:);
                % Get coordintes for relevant atoms
                for jj = 1:size(connected,1);
                    for kk = 1:size(connected,2);
                        coord_find(jj,kk) = {coords_1(find(atom_id_1==connected(jj,kk)),:)};
                    end
                end

                % XYZ format:
                % number of elements
                % comment line
                % <elemente> <X> <Y> <Z>
                % ...

                num_atoms = size(atom_1,1);
                atom_label = char(atom_1);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % NOTE: data normalized to 1 (i.e. relative distance is the same, but overall
                % distance is no longer representative. Must convert back to
                % Angstroms.
                % HOW? ---> use cell volume
                coords_conv = coords_1*cell_param; % Use cell parameter to undo normalization and covert to Angstroms
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                file_xyz = ['Out',num2str(hh),'.xyz']; 

                fileID = fopen(file_xyz, 'wt');
                fprintf(fileID, '%d\n', num_atoms);

                fprintf(fileID, 'NAME\n');

                for ii = 1:size(atom_1,1);
                    fprintf(fileID, '%s\t %i\t %i\t% i\t\n',atom_label(ii,:), coords_conv(ii,:));
                end
                fclose(fileID);

                clear fileID

            end

            %% Convert from xyz to sdf to allow calculation of R3m
            % Run Open Babel to convert
            % see https://openbabel.org/docs/dev/Command-line_tools/babel.html
            % -c = center coordinates
            % -o = specify output format (e.g. -osdf --> output to sdf)

            copyfile(['Out',num2str(hh),'.xyz'], 'Out.xyz'); %Copy file to run in Babel

            %!obabel Out.xyz -c -O Out.sdf
            [yy, zz] = system('obabel Out.xyz -c -O Out.sdf'); % Use OpenBable to convert from .xyz to .sdf
            clear yy zz

            copyfile('Out.sdf', ['Out',num2str(hh),'.sdf']); % Copy output file from Babel to new file name

            %% Calculate R3m
            pathname = cd;
            %file taken in by function

            file_sdf = ['Out',num2str(hh),'.sdf'];
            filename = [pathname,'\', file_sdf];
            % Call other function: R3mCalculate_aut
            [r3m] = R3mCalculate_auto(filename); % Calculate R3m for each molecule

            % r3m_out((hh-minMol+1),aa) = r3m;
            r3m_out(hh,aa) = r3m;

        end

    waitbar(aa/(size(files_in,1)), w, 'Calculating.   Please Wait...'); % update waitbar

    end

end
close(w) % close waitbar 

%% Plot a histogram of calculated R3m values
col = rand(1,3); % Random color
figure;
hold on
histogram(r3m_out, 'FaceColor', col, 'FaceAlpha', 0.5); % Plot histogram
mean_r3m = mean(mean(r3m_out)); % Calculate mean
median_r3m = median(median(r3m_out)); % Calculate median
set(gca, 'FontSize', 12);
% Display Mean on plot
text(0.7,0.95,['Mean R3m = ',num2str(round(mean_r3m,2))], 'Units', 'Normalized','FontSize', 10);
% Display Median on plot
text(0.7,0.9,['Median R3m = ',num2str(round(median_r3m,2))], 'Units', 'Normalized','FontSize', 10);
ylabel('Counts', 'FontSize', 14);
xlabel('R3m', 'FontSize', 14);

print(gcf, 'R3m_histogram', '-djpeg', '-r600') % Save histogram as jpg at 600 dpi

end

