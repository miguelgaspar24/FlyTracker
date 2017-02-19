
% Compute features from raw tracking data.
%
% To compute features, use:
%
%    feat = feat_compute(trk, calib)
%
% where:
%
%    trk         - tracking structure, obtained from previous steps of
%                   tracker (stored in *-track.mat)
%    calib       - calibration structure, obtained from calibration step
%                   (stored in calibration.mat)
%
% returns:
%
%    feat.       - feature structure
%       names    - names of all features
%       units    - units of each feature
%       data     - n_flies x n_frames x n_features matrix containing the 
%                   feature values of all flies at each frame
%
% All features are computed to be independent of resolution (FPS,PPM)
%
function feat = feat_compute(trk, calib)
    % Store video resolution parameters for later normalization.
    pix_per_mm = calib.PPM;
    FPS = calib.FPS;

    % Determine the number of flies to be analyzed.
    n_flies = size(trk.data, 1);
    flies = 1:n_flies;
    comparisons = numel(flies) - 1;
    bud_complete = false(size(flies)); % Flag to avoid repeated comparisons.

    % Names of features to be computed.
    personal_feat = {'vel','ang_vel','min_wing_ang','max_wing_ang',...
                     'mean_wing_length','axis_ratio','fg_body_ratio','contrast'};
    enviro_feat   = {'dist_to_wall'};
    rel_start = numel(personal_feat) + numel(enviro_feat) + 1; % Set the column number to start writing the relative features.
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
                relative_feat(a) = rel_feats(index);
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
                relative_units(a) = rel_units(index);
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
    features = nan(n_flies, n_frames, n_feats);
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
        %%% PERSONAL FEATURES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % VELOCITY & ACCELERATION
        % Smooth out pos data.
        x = track(s,:,1);
        y = track(s,:,2);
        x(2:end-1) = conv(x,smooth_kernel,'valid');
        y(2:end-1) = conv(y,smooth_kernel,'valid');
        % Calculate velocity.
        x_diff = diff(x);
        y_diff = diff(y);
        pos_diff = (x_diff.^2 + y_diff.^2).^0.5;
        vel = [pos_diff(1) pos_diff] * FPS/pix_per_mm;
        features(s,:,1) = vel;
        
        % ANGULAR VELOCITY & ACCELERATION
        % Make angle data continuous.
        ori = track(s,:,3);
        theta2 = ori(2:end);
        theta1 = ori(1:end-1);
        ori_diff = abs(mod(theta1+pi/2-theta2,pi)-pi/2);
        ang_vel = [ori_diff(1) ori_diff] * FPS;
        ang_vel(2:end-1) = conv(ang_vel,smooth_kernel,'valid');
        features(s,:,2) = ang_vel; % should be invariant of turning left or right.

        % WING ANGLES
        ang1 = abs(track(s,:,14));
        ang1(2:end-1) = conv(ang1,smooth_kernel,'valid');
        len1 = track(s,:,15);
        ang2 = abs(track(s,:,16));   
        ang2(2:end-1) = conv(ang2,smooth_kernel,'valid');
        len2 = track(s,:,17);
        mean_wing_lengths = nanmean([len1(:) len2(:)],2);
        min_angles = nanmin([ang1(:) ang2(:)],[],2);
        max_angles = nanmax([ang1(:) ang2(:)],[],2);        
        mean_wing_lengths(2:end-1) = conv(mean_wing_lengths,smooth_kernel,'valid');
        features(s,:,3) = min_angles;
        features(s,:,4) = max_angles;      
        features(s,:,5) = mean_wing_lengths / pix_per_mm;

        % AXIS RATIO
        axis_ratio = track(s,:,4)./track(s,:,5);
        axis_ratio(2:end-1) = conv(axis_ratio,smooth_kernel,'valid');
        features(s,:,6) = axis_ratio;
        
        % FG / BODY RATIO
        fg_body_ratio = track(s,:,7)./track(s,:,6);      
        fg_body_ratio(2:end-1) = conv(fg_body_ratio,smooth_kernel,'valid');
        features(s,:,7) = fg_body_ratio;
        
        % BODY CONTRAST
        body_contrast = track(s,:,8);
        body_contrast(2:end-1) = conv(body_contrast,smooth_kernel,'valid');
        features(s,:,8) = body_contrast;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% RELATIVE TO CHAMBER
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x = min(max(1,round(x)),size(mask,2));
        y = min(max(1,round(y)),size(mask,1));
        inds = sub2ind(size(mask),y,x);
        wall_dist = dists(inds);
        features(s,:,9) = wall_dist / pix_per_mm;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% RELATIVE TO OTHER FLY
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        center = [track(s,:,1)', track(s,:,2)'];
        vec_rot = [cos(ori)', -sin(ori)'];
        buddies = flies(flies~=s); % Determine the buddies for a given fly.
        flag = find(bud_complete);
        % If there is only one fly, then it has no buddies, so we don't
        % need to calculate relative features.
        if buddies == 0
            continue
        end
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
end