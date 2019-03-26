
function writeJAABA(trkname,moviename,trk,feat,calib,n_flies,flies,chamber,s,pair,counts)
    if nargin < 3
        D = load(trkname); trk = D.trk;
    end
    if nargin < 4
        featname = [trkname(1:end-10) '-feat.mat'];
        D = load(featname); feat = D.feat;
    end
    if nargin < 5
        parent_dir = fileparts(fileparths(trkname));
        f_calib = fullfile(parent_dir,'calibration.mat');
        D = load(f_calib); calib = D.calib;
    end

    % Create the experiment directory.
    expdir = [trkname(1:end-10) '_JAABA'];
    if ~exist(expdir,'dir')
      [success1,msg1] = mkdir(expdir);
      if ~success1
        disp(msg1);
        return;
      end
    end
    [~,~,movie_ext] = fileparts(moviename);
    moviefile = fullfile(expdir,['movie' movie_ext]);

    % Write timestamps.
    total_frames = size(trk.data,2);
    n_trkfeat = size(trk.data,3);
    timestamps = (1:total_frames) / calib.FPS;

    % Get wing index.
    wingidx = find(strcmp(trk.names,'wing l x'));

    % Write trx.
    if isfield(trk,'flies_in_chamber')
        if n_flies <= 2
            counter = 1;
            for i=flies
                data = trk.data(i,:,:);
                data = reshape(data,total_frames,n_trkfeat);
                valid_inds = find(~isnan(data(:,1)));
                % Movie in its original folder.
                track.moviename     = moviename;
                % Movie stored in JAABA folder.
                track.moviefile     = moviefile;
                % Video info.
                track.firstframe    = valid_inds(1);
                track.off           = 1-track.firstframe;
                track.endframe      = valid_inds(end);
                track.nframes       = track.endframe - track.firstframe + 1;
                track.fps           = calib.FPS;
                track.pxpermm       = calib.PPM;
                %track.arena %?
                % Fly info.
                track.id            = i;
                track.sex           = 'm';
                % Time stamps.
                frames = track.firstframe:track.endframe;
                track.timestamps    = frames/track.fps;
                track.dt            = diff(track.timestamps);
                % Raw features in pixels.
                track.x             = data(frames,1);
                track.y             = data(frames,2);
                track.theta         = -data(frames,3);
                track.a             = data(frames,4)/4;
                track.b             = data(frames,5)/4;
                track.xwingl        = data(frames,wingidx);
                track.ywingl        = data(frames,wingidx+1);
                track.xwingr        = data(frames,wingidx+2);
                track.ywingr        = data(frames,wingidx+3);
                % Raw features in mm (necessary?).
                track.x_mm          = track.x/track.pxpermm;
                track.y_mm          = track.y/track.pxpermm;
                track.a_mm          = track.a/track.pxpermm;
                track.b_mm          = track.b/track.pxpermm;
                track.theta_mm      = track.theta;

                trx(counter) = track;
                
                counter = counter + 1;
            end

            % Save trx.mat.
            save(fullfile(expdir, ['trx', chamber, '.mat']),'timestamps','trx')

            % Write perframe folder.
            perframedir = fullfile(expdir,['perframe', chamber]);
            if ~exist(perframedir,'dir')
              [success1,msg1] = mkdir(perframedir);
              if ~success1
                disp(msg1);
                return;
              end
            end

        elseif n_flies > 2
            % Write trx.
            counter = 1;
            for i=[s, pair]
                data = trk.data(i,:,:);
                data = reshape(data,total_frames,n_trkfeat);
                valid_inds = find(~isnan(data(:,1)));
                % Movie in its original folder.
                track.moviename     = moviename;
                % Movie stored in JAABA folder.
                track.moviefile     = moviefile;
                % Video info.
                track.firstframe    = valid_inds(1);
                track.off           = 1-track.firstframe;
                track.endframe      = valid_inds(end);
                track.nframes       = track.endframe - track.firstframe + 1;
                track.fps           = calib.FPS;
                track.pxpermm       = calib.PPM;
                % Fly info.
                track.id            = i;
                track.sex           = 'm';
                % Time stamps.
                frames = track.firstframe:track.endframe;
                track.timestamps    = frames/track.fps;
                track.dt            = diff(track.timestamps);
                % Raw features in pixels.
                track.x             = data(frames,1);
                track.y             = data(frames,2);
                track.theta         = -data(frames,3);
                track.a             = data(frames,4)/4;
                track.b             = data(frames,5)/4;
                track.xwingl        = data(frames,wingidx);
                track.ywingl        = data(frames,wingidx+1);
                track.xwingr        = data(frames,wingidx+2);
                track.ywingr        = data(frames,wingidx+3);
                % Raw features in mm (necessary?).
                track.x_mm          = track.x/track.pxpermm;
                track.y_mm          = track.y/track.pxpermm;
                track.a_mm          = track.a/track.pxpermm;
                track.b_mm          = track.b/track.pxpermm;
                track.theta_mm      = track.theta;

                trx(counter) = track;

                counter = counter + 1;
            end

            % Save trx.mat.
            save(fullfile(expdir,['trx', chamber, counts, '.mat']),'timestamps','trx')

            % Write perframe folder.
            perframedir = fullfile(expdir,['perframe', chamber, counts]);
            if ~exist(perframedir,'dir')
              [success1,msg1] = mkdir(perframedir);
              if ~success1
                disp(msg1);
                return;
              end
            end
        end

    % Write trx.
    elseif ~isfield(trk,'flies_in_chamber')
        if n_flies <= 2
            for i=1:n_flies
                data = trk.data(i,:,:);
                data = reshape(data,total_frames,n_trkfeat);
                valid_inds = find(~isnan(data(:,1)));
                % Movie in its original folder.
                track.moviename     = moviename;
                % Movie stored in JAABA folder.
                track.moviefile     = moviefile;
                % Video info.
                track.firstframe    = valid_inds(1);
                track.off           = 1-track.firstframe;
                track.endframe      = valid_inds(end);
                track.nframes       = track.endframe - track.firstframe + 1;
                track.fps           = calib.FPS;
                track.pxpermm       = calib.PPM;
                %track.arena %?
                % Fly info.
                track.id            = i;
                track.sex           = 'm';
                % Time stamps.
                frames = track.firstframe:track.endframe;
                track.timestamps    = frames/track.fps;
                track.dt            = diff(track.timestamps);
                % Raw features in pixels.
                track.x             = data(frames,1);
                track.y             = data(frames,2);
                track.theta         = -data(frames,3);
                track.a             = data(frames,4)/4;
                track.b             = data(frames,5)/4;
                track.xwingl        = data(frames,wingidx);
                track.ywingl        = data(frames,wingidx+1);
                track.xwingr        = data(frames,wingidx+2);
                track.ywingr        = data(frames,wingidx+3);
                % Raw features in mm (necessary?).
                track.x_mm          = track.x/track.pxpermm;
                track.y_mm          = track.y/track.pxpermm;
                track.a_mm          = track.a/track.pxpermm;
                track.b_mm          = track.b/track.pxpermm;
                track.theta_mm      = track.theta;

                trx(i) = track;
                
            end

            % Save trx.mat.
            save(fullfile(expdir, 'trx.mat'),'timestamps','trx')

            % Write perframe folder.
            perframedir = fullfile(expdir,'perframe');
            if ~exist(perframedir,'dir')
              [success1,msg1] = mkdir(perframedir);
              if ~success1
                disp(msg1);
                return;
              end
            end

        elseif n_flies > 2
            % Write trx.
            counter = 1;
            for i=[s, pair]
                data = trk.data(i,:,:);
                data = reshape(data,total_frames,n_trkfeat);
                valid_inds = find(~isnan(data(:,1)));
                % Movie in its original folder.
                track.moviename     = moviename;
                % Movie stored in JAABA folder.
                track.moviefile     = moviefile;
                % Video info.
                track.firstframe    = valid_inds(1);
                track.off           = 1-track.firstframe;
                track.endframe      = valid_inds(end);
                track.nframes       = track.endframe - track.firstframe + 1;
                track.fps           = calib.FPS;
                track.pxpermm       = calib.PPM;
                % Fly info.
                track.id            = i;
                track.sex           = 'm';
                % Time stamps.
                frames = track.firstframe:track.endframe;
                track.timestamps    = frames/track.fps;
                track.dt            = diff(track.timestamps);
                % Raw features in pixels.
                track.x             = data(frames,1);
                track.y             = data(frames,2);
                track.theta         = -data(frames,3);
                track.a             = data(frames,4)/4;
                track.b             = data(frames,5)/4;
                track.xwingl        = data(frames,wingidx);
                track.ywingl        = data(frames,wingidx+1);
                track.xwingr        = data(frames,wingidx+2);
                track.ywingr        = data(frames,wingidx+3);
                % Raw features in mm (necessary?).
                track.x_mm          = track.x/track.pxpermm;
                track.y_mm          = track.y/track.pxpermm;
                track.a_mm          = track.a/track.pxpermm;
                track.b_mm          = track.b/track.pxpermm;
                track.theta_mm      = track.theta;

                trx(counter) = track;

                counter = counter + 1;
            end

            % Save trx.mat.
            save(fullfile(expdir,['trx', counts, '.mat']),'timestamps','trx')

            % Write perframe folder.
            perframedir = fullfile(expdir,['perframe', counts]);
            if ~exist(perframedir,'dir')
              [success1,msg1] = mkdir(perframedir);
              if ~success1
                disp(msg1);
                return;
              end
            end
        end
    end
    
    n_feat = size(feat.data,3);
    for i=1:n_feat
%         data = cell(1,2);
        if n_flies == 1
            data = cell(1,1);
            data = feat.data(:,trx.firstframe:trx.endframe,i);
        else
            data = cell(1,2);
            for s=1:2
                data{s} = feat.data(s,trx(s).firstframe:trx(s).endframe,i);
            end
        end
        units.numerator = cell(1,0);
        units.denominator = cell(1,0);
        save(fullfile(perframedir,[feat.names{i}]),'data','units')
    end

    % Copy/soft-link movie.
    dosoftlink = 1;
    inmoviefile = moviename;
    if dosoftlink
      if exist(moviefile,'file')
        delete(moviefile);
      end
      if isunix
        cmd = sprintf('ln -s %s %s',inmoviefile,moviefile);
        unix(cmd);
        % Test to make sure it worked.
        [status,result] = unix(sprintf('readlink %s',moviefile));
        result = strtrim(result);
        if status ~= 0 || ~strcmp(result,inmoviefile)
          dosoftlink = false;
        end
      elseif ispc
        if exist([moviefile,'.lnk'],'file')
          delete([moviefile,'.lnk']);
        end
        cmd = sprintf('mkshortcut.vbs /target:"%s" /shortcut:"%s"',inmoviefile,moviefile);
        fprintf('Making a Windows shortcut file at "%s" with target "%s"\n',inmoviefile,moviefile);
        system(cmd);
        % Test to make sure that worked.
        if ~exist(moviefile,'file')
            % Try a different softlink method.
            cmd = sprintf('mklink %s %s',moviefile,inmoviefile);
            system(cmd);
            % Test to make sure that worked.
            if ~exist(moviefile,'file')
                dosoftlink = false;
            end
        end
      else
        dosoftlink = false;
      end  
    end

    if ~dosoftlink
      if ispc
        if exist([moviefile,'.lnk'],'file')
          delete([moviefile,'.lnk']);
        end
      end
      if exist(moviefile,'file')
        delete(moviefile);
      end
      disp('WARNING: Could not write softlink');
      %copyfile(inmoviefile,moviefile);
    end

end
