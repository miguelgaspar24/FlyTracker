
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
function feat = feat_compute(f_vid, f_res, options, featfile, trk, calib, recompute)
    % Store video resolution parameters for later normalization.
    pix_per_mm = calib.PPM;
    FPS = calib.FPS;
    
    % Names of features to be computed.
    personal_feat = {'vel','ang_vel','min_wing_ang','max_wing_ang',...
                     'mean_wing_length','axis_ratio','fg_body_ratio','contrast'};
    enviro_feat   = {'dist_to_wall'};
    relative_feat = {'dist_to_other','angle_between','facing_angle','leg_dist'};         
    
    % Units of features to be computed.
    personal_units = {'mm/s','rad/s','rad','rad','mm','ratio','ratio',''};
    enviro_units = {'mm'};
    relative_units = {'mm','rad','rad','mm'};
    
    % Kernel for smoothing output. 
    smooth_kernel = [1 2 1]/4;

    % Set some working variables.
    n_frames = size(trk.data,2);
    n_trkfeat = size(trk.data,3);
    n_feats = numel(personal_feat) + numel(enviro_feat) + numel(relative_feat);
    track = trk.data(:,:,:);
    features = nan(2,n_frames,n_feats);
    
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
    
    
    if isfield(trk,'flies_in_chamber')
        for c=1:numel(trk.flies_in_chamber)
            n_flies = numel(trk.flies_in_chamber{c});
            flies = trk.flies_in_chamber{c};
            bud_complete = zeros(size(flies));
            chamber = sprintf('_c%d', c);
            f = 1;

            if n_flies == 1
                features = nan(1,n_frames,n_feats);
            end
            
            % Compute features for all flies.
            if n_flies <= 2
                for s=flies
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% PERSONAL FEATURES
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % VELOCITY & ACCELERATION
                    %  Smooth out pos data.
                    x = track(s,:,1);
                    y = track(s,:,2);
                    x(2:end-1) = conv(x,smooth_kernel,'valid');
                    y(2:end-1) = conv(y,smooth_kernel,'valid');

                    % Calculate velocity.
                    x_diff = diff(x);
                    y_diff = diff(y);
                    pos_diff = (x_diff.^2 + y_diff.^2).^0.5;
                    vel = [pos_diff(1) pos_diff] * FPS/pix_per_mm;

                    features(f,:,1) = vel;

                    % ANGULAR VELOCITY & ACCELERATION
                    %  Make angle data continuous.
                    ori = track(s,:,3);
                    theta2 = ori(2:end);
                    theta1 = ori(1:end-1);
                    ori_diff = abs(mod(theta1+pi/2-theta2,pi)-pi/2);
                    ang_vel = [ori_diff(1) ori_diff] * FPS;
                    ang_vel(2:end-1) = conv(ang_vel,smooth_kernel,'valid');

                    features(f,:,2) = ang_vel; % Should be invariant of turning left or right.

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

                    features(f,:,3) = min_angles;
                    features(f,:,4) = max_angles;      
                    features(f,:,5) = mean_wing_lengths / pix_per_mm;

                    % AXIS RATIO
                    axis_ratio = track(s,:,4)./track(s,:,5);
                    axis_ratio(2:end-1) = conv(axis_ratio,smooth_kernel,'valid');

                    features(f,:,6) = axis_ratio;

                    % FG / BODY RATIO
                    fg_body_ratio = track(s,:,7)./track(s,:,6);      
                    fg_body_ratio(2:end-1) = conv(fg_body_ratio,smooth_kernel,'valid');

                    features(f,:,7) = fg_body_ratio;

                    % BODY CONTRAST
                    body_contrast = track(s,:,8);
                    body_contrast(2:end-1) = conv(body_contrast,smooth_kernel,'valid');

                    features(f,:,8) = body_contrast;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% RELATIVE TO CHAMBER
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    x = min(max(1,round(x)),size(mask,2));
                    y = min(max(1,round(y)),size(mask,1));
                    inds = sub2ind(size(mask),y,x);
                    wall_dist = dists(inds);

                    features(f,:,9) = wall_dist / pix_per_mm;        

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% RELATIVE TO OTHER FLY
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    pair = flies(flies~=s); % Determine the buddies for a given fly.
                    
                    if size(pair) == 0
                        f = f + 1;
                        continue
                    end
                        
                    center = [track(s,:,1)' track(s,:,2)'];
                    vec_rot = [cos(ori)' -sin(ori)'];

                    %  DISTANCE TO OTHER FLY (+ VELOCITY & ACCELERATION)
                    center2 = [track(pair,:,1)' track(pair,:,2)'];
                    vec_between = center2-center;
                    norm_between = (vec_between(:,1).^2 + vec_between(:,2).^2).^.5;
                    dist_between = norm_between;
                    dist_between(2:end-1) = conv(dist_between,smooth_kernel,'valid');

                    features(1,:,10) = dist_between / pix_per_mm;
                    features(2,:,10) = features(1,:,10);

                    % Compute angles.
                    ori2 = track(pair,:,3);
                    vec_rot2 = [cos(ori2)' -sin(ori2)'];
                    angle_between = acos(dot(vec_rot,vec_rot2,2));
                    vec_between = vec_between./repmat(norm_between,1,2);
                    facing_angle1 = acos(dot(vec_rot,vec_between,2));
                    facing_angle2 = acos(dot(vec_rot2,-vec_between,2));          

                    % ANGLE BETWEEN FLIES
                    angle_between(2:end-1) = conv(angle_between,smooth_kernel,'valid');

                    features(1,:,11) = angle_between;
                    features(2,:,11) = features(1,:,11);

                    % FACING ANGLE 
                    facing_angle1(2:end-1) = conv(facing_angle1,smooth_kernel,'valid');
                    facing_angle2(2:end-1) = conv(facing_angle2,smooth_kernel,'valid');

                    features(1,:,12) = facing_angle1;
                    features(2,:,12) = facing_angle2;

                    % LEG DIST TO OTHER FLY & BODY WING DIST TO OTHER FLY
                    fg_dist = track(s,:,9);
                    fg_dist(2:end-1) = conv(fg_dist,smooth_kernel,'valid');

                    features(1,:,13) = fg_dist / pix_per_mm;
                    features(2,:,13) = features(1,:,13);

                    f = f + 1;
                end
                
                % Store variables in feat structure.
                names = [personal_feat, enviro_feat, relative_feat];
                units = [personal_units, enviro_units, relative_units];

                feat.names = names;
                feat.units = units;
                feat.data = features;

                append = sprintf('-feat_c%d.mat', c);
                savename = [f_res(1:end-10), append];

                % Save xls files.
                if options.save_xls
                    xlsfile = [f_res(1:end-10), '-trackfeat'];
                    if ~exist(xlsfile,'dir') || ~exist([xlsfile '.xlsx'],'file') || recompute
                       if ~exist('trk','var')
                          trk = load(f_res); trk = trk.trk;
                       end
                       if ~exist('feat','var')
                          feat = load(featfile); feat = feat.feat;
                       end
                       names = [trk.names, feat.names];
                       data = nan(n_flies,size(trk.data,2),numel(names));
                       data(:,:,1:size(trk.data,3)) = trk.data(trk.flies_in_chamber{c},:,:);
                       data(:,:,size(trk.data,3)+(1:size(feat.data,3))) = feat.data;
                       writeXls(xlsfile, data, names, trk, n_flies, flies, chamber);
                    end
                end

                % Write JAABA folders.
                if options.save_JAABA
                  JAABA_dir = [f_res(1:end-10) '-JAABA'];
                  if (~exist(JAABA_dir,'dir') || recompute)
                     if ~exist('trk','var')
                         trk = load(f_res);
                         trk = trk.trk;
                     end
                     if ~exist('feat','var')
                         feat = load(featfile);
                         feat = feat.feat;
                     end
                     % Augment features (with log, norms, and derivatives).
                     feat = feat_augment(feat);
                     writeJAABA(f_res, f_vid, trk, feat, calib, n_flies, flies, chamber, s, pair);
                  end
                end

                save(savename, 'feat')

            elseif n_flies > 2
                for s=flies
                    buddies = flies(flies~=s); % Determine the buddies for a given fly.
                    flag = find(bud_complete);
                    for pair=buddies
                        if ismember(pair, flag)
                            continue
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%% PERSONAL FEATURES
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % VELOCITY & ACCELERATION
                        % Smooth out pos data.
                        x = track(s,:,1);
                        y = track(s,:,2);
                        x(2:end-1) = conv(x,smooth_kernel,'valid');
                        y(2:end-1) = conv(y,smooth_kernel,'valid');

                        x_pair = track(pair,:,1);
                        y_pair = track(pair,:,2);
                        x_pair(2:end-1) = conv(x_pair,smooth_kernel,'valid');
                        y_pair(2:end-1) = conv(y_pair,smooth_kernel,'valid');

                        % Calculate velocity.
                        x_diff = diff(x);
                        y_diff = diff(y);
                        pos_diff = (x_diff.^2 + y_diff.^2).^0.5;
                        vel = [pos_diff(1) pos_diff] * FPS/pix_per_mm;

                        x_diff_pair = diff(x_pair);
                        y_diff_pair = diff(y_pair);
                        pos_diff_pair = (x_diff_pair.^2 + y_diff_pair.^2).^0.5;
                        vel_pair = [pos_diff_pair(1) pos_diff_pair] * FPS/pix_per_mm;

                        features(1,:,1) = vel;
                        features(2,:,1) = vel_pair;

                        % ANGULAR VELOCITY & ACCELERATION
                        % Make angle data continuous.
                        ori = track(s,:,3);
                        theta2 = ori(2:end);
                        theta1 = ori(1:end-1);
                        ori_diff = abs(mod(theta1+pi/2-theta2,pi)-pi/2);
                        ang_vel = [ori_diff(1) ori_diff] * FPS;
                        ang_vel(2:end-1) = conv(ang_vel,smooth_kernel,'valid');

                        ori_pair = track(pair,:,3);
                        theta2_pair = ori_pair(2:end);
                        theta1_pair = ori_pair(1:end-1);
                        ori_diff_pair = abs(mod(theta1_pair+pi/2-theta2_pair,pi)-pi/2);
                        ang_vel_pair = [ori_diff_pair(1) ori_diff_pair] * FPS;
                        ang_vel_pair(2:end-1) = conv(ang_vel_pair,smooth_kernel,'valid');

                        features(1,:,2) = ang_vel; % Should be invariant of turning left or right.
                        features(2,:,2) = ang_vel_pair; % Should be invariant of turning left or right.

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

                        ang1_pair = abs(track(pair,:,14));
                        ang1_pair(2:end-1) = conv(ang1_pair,smooth_kernel,'valid');
                        len1_pair = track(pair,:,15);
                        ang2_pair = abs(track(pair,:,16));   
                        ang2_pair(2:end-1) = conv(ang2_pair,smooth_kernel,'valid');
                        len2_pair = track(pair,:,17);
                        mean_wing_lengths_pair = nanmean([len1_pair(:) len2_pair(:)],2);
                        min_angles_pair = nanmin([ang1_pair(:) ang2_pair(:)],[],2);
                        max_angles_pair = nanmax([ang1_pair(:) ang2_pair(:)],[],2);        
                        mean_wing_lengths_pair(2:end-1) = conv(mean_wing_lengths_pair,smooth_kernel,'valid');

                        features(1,:,3) = min_angles;
                        features(1,:,4) = max_angles;      
                        features(1,:,5) = mean_wing_lengths / pix_per_mm;
                        features(2,:,3) = min_angles_pair;
                        features(2,:,4) = max_angles_pair;      
                        features(2,:,5) = mean_wing_lengths_pair / pix_per_mm;

                        % AXIS RATIO
                        axis_ratio = track(s,:,4)./track(s,:,5);
                        axis_ratio(2:end-1) = conv(axis_ratio,smooth_kernel,'valid');

                        axis_ratio_pair = track(pair,:,4)./track(pair,:,5);
                        axis_ratio_pair(2:end-1) = conv(axis_ratio_pair,smooth_kernel,'valid');

                        features(1,:,6) = axis_ratio;
                        features(2,:,6) = axis_ratio_pair;

                        % FG / BODY RATIO
                        fg_body_ratio = track(s,:,7)./track(s,:,6);      
                        fg_body_ratio(2:end-1) = conv(fg_body_ratio,smooth_kernel,'valid');

                        fg_body_ratio_pair = track(pair,:,7)./track(pair,:,6);      
                        fg_body_ratio_pair(2:end-1) = conv(fg_body_ratio_pair,smooth_kernel,'valid');

                        features(1,:,7) = fg_body_ratio;
                        features(2,:,7) = fg_body_ratio_pair;

                        % BODY CONTRAST
                        body_contrast = track(s,:,8);
                        body_contrast(2:end-1) = conv(body_contrast,smooth_kernel,'valid');

                        body_contrast_pair = track(pair,:,8);
                        body_contrast_pair(2:end-1) = conv(body_contrast_pair,smooth_kernel,'valid');

                        features(1,:,8) = body_contrast;
                        features(2,:,8) = body_contrast_pair;

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%% RELATIVE TO CHAMBER
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        x = min(max(1,round(x)),size(mask,2));
                        y = min(max(1,round(y)),size(mask,1));
                        inds = sub2ind(size(mask),y,x);
                        wall_dist = dists(inds);

                        x_pair = min(max(1,round(x_pair)),size(mask,2));
                        y_pair = min(max(1,round(y_pair)),size(mask,1));
                        inds_pair = sub2ind(size(mask),y_pair,x_pair);
                        wall_dist_pair = dists(inds_pair);

                        features(1,:,9) = wall_dist / pix_per_mm;
                        features(2,:,9) = wall_dist_pair / pix_per_mm;

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%% RELATIVE TO OTHER FLY
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % DISTANCE TO OTHER FLY (+ VELOCITY & ACCELERATION)
                        center = [track(s,:,1)', track(s,:,2)'];
                        vec_rot = [cos(ori)', -sin(ori)'];

                        center2 = [track(pair,:,1)', track(pair,:,2)'];
                        vec_between = center2 - center;
                        norm_between = (vec_between(:,1).^2 + vec_between(:,2).^2).^.5;
                        dist_between = norm_between;
                        dist_between(2:end-1) = conv(dist_between,smooth_kernel,'valid');

                        features(1,:,10) = dist_between / pix_per_mm;
                        features(2,:,10) = features(1,:,10);

                        % Compute angles.
                        ori2 = track(pair,:,3);
                        vec_rot2 = [cos(ori2)' -sin(ori2)'];
                        angle_between = acos(dot(vec_rot,vec_rot2,2));
                        vec_between = vec_between./repmat(norm_between,1,2);
                        facing_angle1 = acos(dot(vec_rot,vec_between,2));
                        facing_angle2 = acos(dot(vec_rot2,-vec_between,2));          

                        % ANGLE BETWEEN FLIES
                        angle_between(2:end-1) = conv(angle_between,smooth_kernel,'valid');

                        features(1,:,11) = angle_between;
                        features(2,:,11) = features(1,:,11);

                        % FACING ANGLE
                        facing_angle1(2:end-1) = conv(facing_angle1,smooth_kernel,'valid');
                        facing_angle2(2:end-1) = conv(facing_angle2,smooth_kernel,'valid');

                        features(1,:,12) = facing_angle1;
                        features(2,:,12) = facing_angle2;

                        % LEG DIST TO OTHER FLY & BODY WING DIST TO OTHER FLY
                        fg_dist = track(s,:,9);
                        fg_dist(2:end-1) = conv(fg_dist,smooth_kernel,'valid');

                        features(1,:,13) = fg_dist / pix_per_mm;
                        features(2,:,13) = features(1,:,13);

                        % Store variables in feat structure.
                        feat.data = features;

                        names = [personal_feat, enviro_feat, relative_feat];
                        units = [personal_units, enviro_units, relative_units];

                        feat.names = names;
                        feat.units = units;
                        
                        append = sprintf('-feat_c%d_%d_%d.mat', c, s, pair);
                        savename = [f_res(1:end-10), append];
                        if s < 10 && pair < 10
                            counts = savename(end-7:end-4);
                        elseif s >= 10 && pair >= 10
                            counts = savename(end-9:end-4);
                        else
                            counts = savename(end-8:end-4);
                        end

                        % Save xls files.
                        if options.save_xls
                            xlsfile = [f_res(1:end-10), '-trackfeat'];
                            if ~exist(xlsfile,'dir') || ~exist([xlsfile '.xlsx'],'file') || recompute
                               if ~exist('trk','var')
                                  trk = load(f_res); trk = trk.trk;
                               end
                               if ~exist('feat','var')
                                  feat = load(featfile); feat = feat.feat;
                               end
                               names = [trk.names, feat.names];
                               data = nan(2,size(trk.data,2),numel(names));
                               data(:,:,1:size(trk.data,3)) = trk.data([s, pair], :, :);
                               data(:,:,size(trk.data,3)+(1:size(feat.data,3))) = feat.data;
                               writeXls(xlsfile, data, names, trk, n_flies, [], chamber, s, pair, counts);
                            end
                        end

                        % Write JAABA folders.
                        if options.save_JAABA
                          JAABA_dir = [f_res(1:end-10) '-JAABA'];
                          if (~exist(JAABA_dir,'dir') || recompute)
                             if ~exist('trk','var')
                                 trk = load(f_res); trk = trk.trk;
                             end
                             if ~exist('feat','var')
                                 feat = load(featfile); feat = feat.feat;
                             end
                             % Augment features (with log, norms, and derivatives).
                             feat = feat_augment(feat);
                             writeJAABA(f_res, f_vid, trk, feat, calib, n_flies, [], chamber, s, pair, counts);
                          end
                        end

                        save(savename, 'feat')
                    end

                    % Update the flag.
                    bud_complete(s) = 1;
                end
            end
        end


    elseif ~isfield(trk,'flies_in_chamber')
        n_flies = size(trk.data, 1);
        flies = 1:n_flies;
        bud_complete = zeros(size(flies));
        
        if n_flies == 1
            features = nan(1,n_frames,n_feats);
        end
        
        % Compute features for all flies.
        if n_flies <= 2
            for s=flies
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% PERSONAL FEATURES
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % VELOCITY & ACCELERATION
                %  Smooth out pos data.
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
                %  Make angle data continuous.
                ori = track(s,:,3);
                theta2 = ori(2:end);
                theta1 = ori(1:end-1);
                ori_diff = abs(mod(theta1+pi/2-theta2,pi)-pi/2);
                ang_vel = [ori_diff(1) ori_diff] * FPS;
                ang_vel(2:end-1) = conv(ang_vel,smooth_kernel,'valid');

                features(s,:,2) = ang_vel; % Should be invariant of turning left or right.

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
                bud_complete = zeros(size(flies));
                buddies = flies(flies~=s); % Determine the buddies for a given fly.
                flag = find(bud_complete);
                for pair=buddies
                    if ismember(pair, flag) % Avoid repeating calculations for a pair that has already been computed.
                        continue
                    end
                    
                    if size(buddies) == 0 % If there is only one fly, then there are no buddies, and therefore no relative features to compute.
                        continue
                    end

                    center = [track(s,:,1)' track(s,:,2)'];
                    vec_rot = [cos(ori)' -sin(ori)'];

                    %  DISTANCE TO OTHER FLY (+ VELOCITY & ACCELERATION)
                    center2 = [track(pair,:,1)' track(pair,:,2)'];
                    vec_between = center2-center;
                    norm_between = (vec_between(:,1).^2 + vec_between(:,2).^2).^.5;
                    dist_between = norm_between;
                    dist_between(2:end-1) = conv(dist_between,smooth_kernel,'valid');

                    features(1,:,10) = dist_between / pix_per_mm;
                    features(2,:,10) = features(1,:,10);

                    % Compute angles.
                    ori2 = track(pair,:,3);
                    vec_rot2 = [cos(ori2)' -sin(ori2)'];
                    angle_between = acos(dot(vec_rot,vec_rot2,2));
                    vec_between = vec_between./repmat(norm_between,1,2);
                    facing_angle1 = acos(dot(vec_rot,vec_between,2));
                    facing_angle2 = acos(dot(vec_rot2,-vec_between,2));          

                    % ANGLE BETWEEN FLIES
                    angle_between(2:end-1) = conv(angle_between,smooth_kernel,'valid');

                    features(1,:,11) = angle_between;
                    features(2,:,11) = features(1,:,11);

                    % FACING ANGLE 
                    facing_angle1(2:end-1) = conv(facing_angle1,smooth_kernel,'valid');
                    facing_angle2(2:end-1) = conv(facing_angle2,smooth_kernel,'valid');

                    features(1,:,12) = facing_angle1;
                    features(2,:,12) = facing_angle2;

                    % LEG DIST TO OTHER FLY & BODY WING DIST TO OTHER FLY
                    fg_dist = track(s,:,9);
                    fg_dist(2:end-1) = conv(fg_dist,smooth_kernel,'valid');

                    features(1,:,13) = fg_dist / pix_per_mm;
                    features(2,:,13) = features(1,:,13);
                    
                    bud_complete(s) = 1;
                end
            end
            
            % Store variables in feat structure.
            names = [personal_feat, enviro_feat, relative_feat];
            units = [personal_units, enviro_units, relative_units];

            feat.names = names;
            feat.units = units;
            feat.data = features;

            savename = [f_res(1:end-10), '-feat.mat'];

             % Save xls files.
            if options.save_xls
                xlsfile = [f_res(1:end-10), '-trackfeat'];
                if ~exist(xlsfile,'dir') || ~exist([xlsfile '.xlsx'],'file') || recompute
                   if ~exist('trk','var')
                      trk = load(f_res); trk = trk.trk;
                   end
                   if ~exist('feat','var')
                      feat = load(featfile); feat = feat.feat;
                   end
                   names = [trk.names, feat.names];
                   data = nan(n_flies,size(trk.data,2),numel(names));
                   data(:,:,1:size(trk.data,3)) = trk.data;
                   data(:,:,size(trk.data,3)+(1:size(feat.data,3))) = feat.data;
                   writeXls(xlsfile, data, names, trk, n_flies);
                end
            end

            % Write JAABA folders.
            if options.save_JAABA
              JAABA_dir = [f_res(1:end-10) '-JAABA'];
              if (~exist(JAABA_dir,'dir') || recompute)
                 if ~exist('trk','var')
                     trk = load(f_res);
                     trk = trk.trk;
                 end
                 if ~exist('feat','var')
                     feat = load(featfile);
                     feat = feat.feat;
                 end
                 % Augment features (with log, norms, and derivatives).
                 feat = feat_augment(feat);
                 writeJAABA(f_res, f_vid, trk, feat, calib, n_flies);
              end
            end

            save(savename, 'feat')
            
        elseif n_flies > 2
            for s=1:n_flies
                buddies = flies(flies~=s); % Determine the buddies for a given fly.
                flag = find(bud_complete);
                for pair=buddies
                    if ismember(pair, flag)
                        continue
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% PERSONAL FEATURES
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % VELOCITY & ACCELERATION
                    % Smooth out pos data.
                    x = track(s,:,1);
                    y = track(s,:,2);
                    x(2:end-1) = conv(x,smooth_kernel,'valid');
                    y(2:end-1) = conv(y,smooth_kernel,'valid');

                    x_pair = track(pair,:,1);
                    y_pair = track(pair,:,2);
                    x_pair(2:end-1) = conv(x_pair,smooth_kernel,'valid');
                    y_pair(2:end-1) = conv(y_pair,smooth_kernel,'valid');

                    % Calculate velocity.
                    x_diff = diff(x);
                    y_diff = diff(y);
                    pos_diff = (x_diff.^2 + y_diff.^2).^0.5;
                    vel = [pos_diff(1) pos_diff] * FPS/pix_per_mm;

                    x_diff_pair = diff(x_pair);
                    y_diff_pair = diff(y_pair);
                    pos_diff_pair = (x_diff_pair.^2 + y_diff_pair.^2).^0.5;
                    vel_pair = [pos_diff_pair(1) pos_diff_pair] * FPS/pix_per_mm;

                    features(1,:,1) = vel;
                    features(2,:,1) = vel_pair;

                    % ANGULAR VELOCITY & ACCELERATION
                    % Make angle data continuous.
                    ori = track(s,:,3);
                    theta2 = ori(2:end);
                    theta1 = ori(1:end-1);
                    ori_diff = abs(mod(theta1+pi/2-theta2,pi)-pi/2);
                    ang_vel = [ori_diff(1) ori_diff] * FPS;
                    ang_vel(2:end-1) = conv(ang_vel,smooth_kernel,'valid');

                    ori_pair = track(pair,:,3);
                    theta2_pair = ori_pair(2:end);
                    theta1_pair = ori_pair(1:end-1);
                    ori_diff_pair = abs(mod(theta1_pair+pi/2-theta2_pair,pi)-pi/2);
                    ang_vel_pair = [ori_diff_pair(1) ori_diff_pair] * FPS;
                    ang_vel_pair(2:end-1) = conv(ang_vel_pair,smooth_kernel,'valid');

                    features(1,:,2) = ang_vel; % Should be invariant of turning left or right.
                    features(2,:,2) = ang_vel_pair; % Should be invariant of turning left or right.

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

                    ang1_pair = abs(track(pair,:,14));
                    ang1_pair(2:end-1) = conv(ang1_pair,smooth_kernel,'valid');
                    len1_pair = track(pair,:,15);
                    ang2_pair = abs(track(pair,:,16));   
                    ang2_pair(2:end-1) = conv(ang2_pair,smooth_kernel,'valid');
                    len2_pair = track(pair,:,17);
                    mean_wing_lengths_pair = nanmean([len1_pair(:) len2_pair(:)],2);
                    min_angles_pair = nanmin([ang1_pair(:) ang2_pair(:)],[],2);
                    max_angles_pair = nanmax([ang1_pair(:) ang2_pair(:)],[],2);        
                    mean_wing_lengths_pair(2:end-1) = conv(mean_wing_lengths_pair,smooth_kernel,'valid');

                    features(1,:,3) = min_angles;
                    features(1,:,4) = max_angles;      
                    features(1,:,5) = mean_wing_lengths / pix_per_mm;
                    features(2,:,3) = min_angles_pair;
                    features(2,:,4) = max_angles_pair;      
                    features(2,:,5) = mean_wing_lengths_pair / pix_per_mm;

                    % AXIS RATIO
                    axis_ratio = track(s,:,4)./track(s,:,5);
                    axis_ratio(2:end-1) = conv(axis_ratio,smooth_kernel,'valid');

                    axis_ratio_pair = track(pair,:,4)./track(pair,:,5);
                    axis_ratio_pair(2:end-1) = conv(axis_ratio_pair,smooth_kernel,'valid');

                    features(1,:,6) = axis_ratio;
                    features(2,:,6) = axis_ratio_pair;

                    % FG / BODY RATIO
                    fg_body_ratio = track(s,:,7)./track(s,:,6);      
                    fg_body_ratio(2:end-1) = conv(fg_body_ratio,smooth_kernel,'valid');

                    fg_body_ratio_pair = track(pair,:,7)./track(pair,:,6);      
                    fg_body_ratio_pair(2:end-1) = conv(fg_body_ratio_pair,smooth_kernel,'valid');

                    features(1,:,7) = fg_body_ratio;
                    features(2,:,7) = fg_body_ratio_pair;

                    % BODY CONTRAST
                    body_contrast = track(s,:,8);
                    body_contrast(2:end-1) = conv(body_contrast,smooth_kernel,'valid');

                    body_contrast_pair = track(pair,:,8);
                    body_contrast_pair(2:end-1) = conv(body_contrast_pair,smooth_kernel,'valid');

                    features(1,:,8) = body_contrast;
                    features(2,:,8) = body_contrast_pair;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% RELATIVE TO CHAMBER
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    x = min(max(1,round(x)),size(mask,2));
                    y = min(max(1,round(y)),size(mask,1));
                    inds = sub2ind(size(mask),y,x);
                    wall_dist = dists(inds);

                    x_pair = min(max(1,round(x_pair)),size(mask,2));
                    y_pair = min(max(1,round(y_pair)),size(mask,1));
                    inds_pair = sub2ind(size(mask),y_pair,x_pair);
                    wall_dist_pair = dists(inds_pair);

                    features(1,:,9) = wall_dist / pix_per_mm;
                    features(2,:,9) = wall_dist_pair / pix_per_mm;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% RELATIVE TO OTHER FLY
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % DISTANCE TO OTHER FLY (+ VELOCITY & ACCELERATION)
                    center = [track(s,:,1)', track(s,:,2)'];
                    vec_rot = [cos(ori)', -sin(ori)'];

                    center2 = [track(pair,:,1)', track(pair,:,2)'];
                    vec_between = center2 - center;
                    norm_between = (vec_between(:,1).^2 + vec_between(:,2).^2).^.5;
                    dist_between = norm_between;
                    dist_between(2:end-1) = conv(dist_between,smooth_kernel,'valid');

                    features(1,:,10) = dist_between / pix_per_mm;
                    features(2,:,10) = features(1,:,10);

                    % Compute angles.
                    ori2 = track(pair,:,3);
                    vec_rot2 = [cos(ori2)' -sin(ori2)'];
                    angle_between = acos(dot(vec_rot,vec_rot2,2));
                    vec_between = vec_between./repmat(norm_between,1,2);
                    facing_angle1 = acos(dot(vec_rot,vec_between,2));
                    facing_angle2 = acos(dot(vec_rot2,-vec_between,2));          

                    % ANGLE BETWEEN FLIES
                    angle_between(2:end-1) = conv(angle_between,smooth_kernel,'valid');

                    features(1,:,11) = angle_between;
                    features(2,:,11) = features(1,:,11);

                    % FACING ANGLE
                    facing_angle1(2:end-1) = conv(facing_angle1,smooth_kernel,'valid');
                    facing_angle2(2:end-1) = conv(facing_angle2,smooth_kernel,'valid');

                    features(1,:,12) = facing_angle1;
                    features(2,:,12) = facing_angle2;

                    % LEG DIST TO OTHER FLY & BODY WING DIST TO OTHER FLY
                    fg_dist = track(s,:,9);
                    fg_dist(2:end-1) = conv(fg_dist,smooth_kernel,'valid');

                    features(1,:,13) = fg_dist / pix_per_mm;
                    features(2,:,13) = features(1,:,13);

                    % Store variables in feat structure.
                    feat.data = features;

                    names = [personal_feat, enviro_feat, relative_feat];
                    units = [personal_units, enviro_units, relative_units];

                    feat.names = names;
                    feat.units = units;

                    append = sprintf('-feat_%d_%d.mat', s, pair);
                    savename = [f_res(1:end-10), append];
                    counts = savename(end-7:end-4);

                    % Save xls files.
                    if options.save_xls
                        xlsfile = [f_res(1:end-10), '-trackfeat'];
                        if ~exist(xlsfile,'dir') || ~exist([xlsfile '.xlsx'],'file') || recompute
                           if ~exist('trk','var')
                              trk = load(f_res); trk = trk.trk;
                           end
                           if ~exist('feat','var')
                              feat = load(featfile); feat = feat.feat;
                           end
                           names = [trk.names, feat.names];
                           data = nan(2,size(trk.data,2),numel(names));
                           data(:,:,1:size(trk.data,3)) = trk.data([s, pair], :, :);
                           data(:,:,size(trk.data,3)+(1:size(feat.data,3))) = feat.data;
                           writeXls(xlsfile, data, names, trk, n_flies, [], [], s, pair, counts);
                        end
                    end

                    % Write JAABA folders.
                    if options.save_JAABA
                      JAABA_dir = [f_res(1:end-10) '-JAABA'];
                      if (~exist(JAABA_dir,'dir') || recompute)
                         if ~exist('trk','var')
                             trk = load(f_res); trk = trk.trk;
                         end
                         if ~exist('feat','var')
                             feat = load(featfile); feat = feat.feat;
                         end
                         % Augment features (with log, norms, and derivatives).
                         feat = feat_augment(feat);
                         writeJAABA(f_res, f_vid, trk, feat, calib, n_flies, [], [], s, pair, counts);
                      end
                    end

                    save(savename, 'feat')
                end

                % Update the flag.
                bud_complete(s) = 1;
            end
        end
    end
end
