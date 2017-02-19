
%
% Adapted from: feat_compute.m (Eyrún's original version of the script).
%
% This is the main file for the FlyTracker patch. This is the only file
% that needs to be open in MATLAB, and by running it, all companion files
% will be run as well. The only manual step required before running the
% patch is to set the 'folderspath' variable (line 15) to the path where
% all FlyTracker output folders that need to be patched are located on your
% computer.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CHANGE ME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folderspath = 'C:\Users\Miguel\Desktop\Cecilia';

%%%% CHANGE ME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% If we start MATLAB in a different path, change directory to our current
% 'folderspath' and load the calibration file that should also be found there.
cd(folderspath);
load('calibration.mat');

% Find all the files we have in our 'folderspath', keep only the ones that
% are folders, and exclude everything that starts with a dot (in order to
% avoid possible hidden files, as well as '.' and '..' folders).
dirlist = dir(folderspath);
dirlist = dirlist([dirlist.isdir]==1);
dirlist = dirlist(arrayfun(@(x) x.name(1), dirlist) ~= '.');

% Now loop through all the folders that FlyTracker generated, and patch all
% necessary files found in each of those folders.
error_list = {};
for a = 1:length(dirlist)
    try
        % Get information about the location of .mat files, then look for
        % and load both -feat.mat and -track.mat files into the workspace.
        dirinfo = dirlist(a);
        subdir = what([dirinfo.folder, '\', dirinfo.name]);
        matfiles = subdir.mat;
        folder = subdir.path;
        cd(folder);
        [~, terminal_folder, ~] = fileparts(folder);
        feat_finder = ~cellfun('isempty', strfind(matfiles, '-feat.mat'));
        feat_index = find(feat_finder);
        featname = matfiles{feat_index};
        trk_finder = ~cellfun('isempty', strfind(matfiles, '-track.mat'));
        trk_index = find(trk_finder);
        trkname = [folder, '\', matfiles{trk_index}];
        moviename = [folder, '\', terminal_folder, '.avi'];
        load(featname);
        load(trkname);

        % Store video resolution parameters for later normalization.
        pix_per_mm = calib.PPM;

        % Determine the number of flies to be analyzed.
        n_flies = size(trk.data, 1);
        flies = 1:n_flies;
        comparisons = numel(flies) - 1;
        bud_complete = false(size(flies)); % Flag to avoid repeated comparisons.

        % Names of features to be computed.
        personal_feat = {'vel','ang_vel','min_wing_ang','max_wing_ang',...
                         'mean_wing_length','axis_ratio','fg_body_ratio','contrast'};
        enviro_feat   = {'dist_to_wall'};
        rel_start = 10; % Set the column number to start writing the relative features.
        rel_feat_number = 4; % Number of relative features we want to calculate.
        rel_feat_length = rel_feat_number * comparisons; % Total number of columns we will need.
        relative_feat = cell(1, rel_feat_length);
        % Loop through the number of comparisons we have to make (dependent on the
        % number of flies), numbering the different features accordingly, storing
        % their names in a variable (rel_feats), and then order them back in
        % relative_feat.
        for comparison = 1:comparisons
            rel_feat1 = ['dist_to_other', int2str(comparison)];
            rel_feat2 = ['angle_between', int2str(comparison)];
            rel_feat3 = ['facing_angle', int2str(comparison)];
            rel_feat4 = ['leg_dist', int2str(comparison)];
            rel_feats = {rel_feat1, rel_feat2, rel_feat3, rel_feat4};
            index = 1;
            if n_flies == 2
                for a = 1:rel_feat_length
                    relative_feat((a-1)+comparison) = rel_feats(index);
                    index = index + 1;
                end
            elseif n_flies > 2
                for a = linspace(0, rel_feat_length-comparisons, 4)
                    relative_feat(a+comparison) = rel_feats(index);
                    index = index + 1;
                end
            end
        end

        % Units of features to be computed.
        personal_units = {'mm/s','rad/s','rad','rad','mm','ratio','ratio',''};
        enviro_units = {'mm'};
        % Perform the same computations as above to determine the number, name, and
        % placement of the relative units.
        relative_units = cell(1, rel_feat_length);

        for comparison = 1:comparisons
            rel_unit1 = 'mm';
            rel_unit2 = 'rad';
            rel_unit3 = 'rad';
            rel_unit4 = 'mm';
            rel_units = {rel_unit1, rel_unit2, rel_unit3, rel_unit4};
            index = 1;
            if n_flies == 2
                for a = 1:rel_feat_length
                    relative_units((a-1)+comparison) = rel_units(index);
                    index = index + 1;
                end
            elseif n_flies > 2
                for a = linspace(0, rel_feat_length-comparisons, 4)
                    relative_units(a+comparison) = rel_units(index);
                    index = index + 1;
                end
            end
        end

        % Kernel for smoothing output.
        smooth_kernel = [1 2 1]/4;
        % Set some working variables.
        n_frames = size(trk.data,2);
        n_trkfeat = size(trk.data,3);
        n_feats = numel(personal_feat) + (n_flies >= 2) * numel(relative_feat) + numel(enviro_feat);
        track = trk.data(:,:,:);
        present_data = feat.data(:,:,:);
        % Avoid to keep growing dataset if we rerun this script on a dataset that
        % has already been patched. This will stop the rest of the script from
        % running.
        if size(present_data, 3) > numel(personal_feat) + numel(enviro_feat)
            error('This dataset has already been patched. Try another!')
        end
        missing_data = nan(size(feat.data,1), size(feat.data,2), numel(relative_feat));
        features = cat(3, present_data, missing_data);
        column_number = 0:(numel(relative_feat)/comparisons)-1;


        % Compute distance to chambers for all pixels.
        mask = zeros(size(calib.mask));
        for i=1:numel(calib.masks)
            mask = mask | calib.masks{i};
        end
        dists = bwdist(1-mask);


        % Interpolate track values.
        for s=1:size(trk.data,1)
            for f_ind=1:n_trkfeat
                vec = track(s,:,f_ind);
                invalid = isnan(vec);
                cc = bwconncomp(invalid);
                for c=1:cc.NumObjects
                    fr_start = cc.PixelIdxList{c}(1)-1;
                    fr_end   = cc.PixelIdxList{c}(end)+1;
                    frs = fr_start:fr_end;
                    if fr_start < 1 || fr_end > n_frames
                        continue
                    end
                    piece = (vec(fr_end)-vec(fr_start))/(numel(frs)-1);
                    coeffs = 0:(numel(frs)-1);
                    vec(frs) = vec(fr_start) + coeffs * piece;
                end
                track(s,:,f_ind) = vec;
            end
        end

        % Compute features for all flies.
        for s=1:n_flies

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% RELATIVE TO CHAMBER
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Smooth out pos data.
            x = track(s,:,1);
            y = track(s,:,2);
            x(2:end-1) = conv(x,smooth_kernel,'valid');
            y(2:end-1) = conv(y,smooth_kernel,'valid');
            x = min(max(1,round(x)),size(mask,2));
            y = min(max(1,round(y)),size(mask,1));
            inds = sub2ind(size(mask),y,x);
            wall_dist = dists(inds);
            features(s,:,9) = wall_dist / pix_per_mm;



            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% RELATIVE TO OTHER FLY
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            center = [track(s,:,1)', track(s,:,2)'];
            ori = track(s,:,3);
            vec_rot = [cos(ori)', -sin(ori)'];
            buddies = flies(flies~=s); % Determine the buddies for a given fly.
            flag = find(bud_complete);
            % If 1-2 has been computed, avoid wasting time computing 2-1, and so on.
            for pair = buddies
                if ismember(pair, flag)
                    continue
                end
                % DISTANCE TO OTHER FLY (+ VELOCITY & ACCELERATION)
                s_disttoother_column = rel_start + comparisons * column_number(1);
                pair_disttoother_column = rel_start + comparisons * column_number(1);
                center2 = [track(pair,:,1)', track(pair,:,2)'];
                vec_between = center2 - center;
                norm_between = (vec_between(:,1).^2 + vec_between(:,2).^2).^.5;
                dist_between = norm_between;
                dist_between(2:end-1) = conv(dist_between,smooth_kernel,'valid');
                % Keep searching for an empty column to write in, in order to avoid
                % overwriting data.
                while ~all(isnan(features(s,:,s_disttoother_column)))
                    s_disttoother_column = s_disttoother_column + 1;
                end
                features(s,:,s_disttoother_column) = dist_between / pix_per_mm;
                while ~all(isnan(features(pair,:,pair_disttoother_column)))
                    pair_disttoother_column = pair_disttoother_column + 1;
                end
                features(pair,:,pair_disttoother_column) = features(s,:,s_disttoother_column);

                % Compute angles.
                ori2 = track(pair,:,3);
                vec_rot2 = [cos(ori2)' -sin(ori2)'];
                angle_between = acos(dot(vec_rot,vec_rot2,2));
                vec_between = vec_between./repmat(norm_between,1,2);
                facing_angle1 = acos(dot(vec_rot,vec_between,2));
                facing_angle2 = acos(dot(vec_rot2,-vec_between,2));          

                % ANGLE BETWEEN FLIES
                s_anglebetween_column = rel_start + comparisons * column_number(2);
                pair_anglebetween_column = rel_start + comparisons * column_number(2);
                angle_between(2:end-1) = conv(angle_between,smooth_kernel,'valid');
                while ~all(isnan(features(s,:,s_anglebetween_column)))
                    s_anglebetween_column = s_anglebetween_column + 1;
                end
                features(s,:,s_anglebetween_column) = angle_between;
                while ~all(isnan(features(pair,:,pair_anglebetween_column)))
                    pair_anglebetween_column = pair_anglebetween_column + 1;
                end
                features(pair,:,pair_anglebetween_column) = features(s,:,s_anglebetween_column);

                % FACING ANGLE
                s_facingangle_column = rel_start + comparisons * column_number(3);
                pair_facingangle_column = rel_start + comparisons * column_number(3);
                facing_angle1(2:end-1) = conv(facing_angle1,smooth_kernel,'valid');
                facing_angle2(2:end-1) = conv(facing_angle2,smooth_kernel,'valid');
                while ~all(isnan(features(s,:,s_facingangle_column)))
                    s_facingangle_column = s_facingangle_column + 1;
                end
                features(s,:,s_facingangle_column) = facing_angle1;
                while ~all(isnan(features(pair,:,pair_facingangle_column)))
                    pair_facingangle_column = pair_facingangle_column + 1;
                end
                features(pair,:,pair_facingangle_column) = facing_angle2;

                % LEG DIST TO OTHER FLY & BODY WING DIST TO OTHER FLY
                s_legdist_column = rel_start + comparisons * column_number(4);
                pair_legdist_column = rel_start + comparisons * column_number(4);
                fg_dist = track(s,:,9);
                fg_dist(2:end-1) = conv(fg_dist,smooth_kernel,'valid');
                while ~all(isnan(features(s,:,s_legdist_column)))
                    s_legdist_column = s_legdist_column + 1;
                end
                features(s,:,s_legdist_column) = fg_dist / pix_per_mm;
                while ~all(isnan(features(pair,:,pair_legdist_column)))
                    pair_legdist_column = pair_legdist_column + 1;
                end
                features(pair,:,pair_legdist_column) = features(s,:,s_legdist_column);

                % Store variables in feat structure.
                feat.data = features;
            end

            % Store variables in feat structure.
            names = [personal_feat, enviro_feat, relative_feat];   
            units = [personal_units, enviro_units, relative_units];

            feat.names = names;
            feat.units = units;

            % Update the flag.
            bud_complete(s) = 1;
        end

        % Save the newly computed relative features.
        save(featname, 'feat')

        % Run the next components of the patching process.
        run('feat_augment_patch.m')
        run('writeJAABA_patch.m')
        run('writeXls_patch.m')

        % Clear 'trk and 'feat' variables to avoid possible conflicts when
        % loading the 'trk' and 'feat' variables of the next folder.
        clearvars('trk');
        disp(' ')
        disp(['trk variable of folder ', terminal_folder, ' cleared'])
        clearvars('feat');
        disp(['feat variable of folder ', terminal_folder, ' cleared'])
        
    % Anytime an error occurs, create a .txt log, and append that error's
    % message to that file. This allows for the patch to run continously
    % through all the folders without being interrupted anytime an error
    % would be triggered. The user can then look at the errors post-run.
    catch err
        fid = fopen('error_log.txt','a+t', 'n');
        fprintf(fid, '%s', err.getReport('extended', 'hyperlinks','off'));
        fprintf(fid, '\n\n\n');
        fclose(fid);
    end
    
    % Get a list of all the .txt error logs generated during patching.
    [~,list] = system('dir /S/B *.txt');
    list = textscan(list, '%s', 'Delimiter', '\n');
    error_list = [error_list; list{:}];
end

% Calculate the number of errors that have occurred.
n_errors = length(error_list);
log = string();

% If there is at least one error in the list, and that element is not a
% File Not Found, then we look for the folders they are in, and display a
% warning message to the user, that reports the number of errors that have
% occurred, and all the folders to look in for the respective error logs.
%error_matrix = cell2mat(error_list);


% Otherwise, display a message box informing the user that the patching
% proccess completed successfully, and without any errors along the way.
if all(strcmp(error_list, 'File Not Found'))
    info = msgbox({'Patch successful. All relevant files restored.';'';...
                   'Number of errors found while patching: 0';''},...
                   'Success!',...
                   'Help',...
                   'help');
               
else
    for error = 1:n_errors
        [~, err_folder, ~] = fileparts(fileparts(error_list{error}));
        log(end+1) = err_folder;
    end

    log = char(log(2:end));
    log_list = [];
    for entry = 1:n_errors
        log_list = [log_list, ',', log(:,:,entry)];
    end

    entry_list = char(strsplit(log_list, ','));
    
    warning = errordlg({'Careful! Looks like some errors occurred during the patching process...';'';...
                        ['Number of errors found while patching: ', num2str(n_errors)];'';...
                        'Please check the following folders for the respective error logs:';entry_list},...
                        'Error Summary',...
                        'modal');
end
