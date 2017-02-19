
% Adapted from: writeJAABA.m (Eyrún's original version of the script).

try
    % Create the experiment directory.
    expdir = [trkname(1:end-10), '_JAABA'];
    if ~exist(expdir,'dir')
      [success1,msg1] = mkdir(expdir);
      if ~success1
        disp(msg1);
        return;
      end
    end
    [~,~,movie_ext] = fileparts(moviename);
    moviefile = fullfile(expdir,['movie', movie_ext]);

    % Write timestamps.
    total_frames = size(trk.data,2);
    n_trkfeat = size(trk.data,3);
    n_flies = size(trk.data,1);
    timestamps = (1:total_frames) / calib.FPS;

    % Get wing index.
    wingidx = find(strcmp(trk.names,'wing l x'));

    % Write trx.
    for i=1:n_flies
        data = trk.data(i,:,:);
        data = reshape(data,total_frames,n_trkfeat);
        valid_inds = find(~isnan(data(:,1)));
        % Movie in its original folder.
        newtrack.moviename     = moviename;
        % Movie stored in JAABA folder.
        newtrack.moviefile     = moviefile;
        % Video info.
        newtrack.firstframe    = valid_inds(1);
        newtrack.off           = 1-newtrack.firstframe;
        newtrack.endframe      = valid_inds(end);
        newtrack.nframes       = newtrack.endframe - newtrack.firstframe + 1;
        newtrack.fps           = calib.FPS;
        newtrack.pxpermm       = calib.PPM;
        % Fly info.
        newtrack.id            = i;
        newtrack.sex           = 'm';
        % Time stamps.
        frames = newtrack.firstframe:newtrack.endframe;
        newtrack.timestamps    = frames/newtrack.fps;
        newtrack.dt            = diff(newtrack.timestamps);
        % Raw features in pixels.
        newtrack.x             = data(frames,1);
        newtrack.y             = data(frames,2);
        newtrack.theta         = -data(frames,3);
        newtrack.a             = data(frames,4)/4;
        newtrack.b             = data(frames,5)/4;
        newtrack.xwingl        = data(frames,wingidx);
        newtrack.ywingl        = data(frames,wingidx+1);
        newtrack.xwingr        = data(frames,wingidx+2);
        newtrack.ywingr        = data(frames,wingidx+3);
        % Raw features in mm. (necessary?)
        newtrack.x_mm          = newtrack.x/newtrack.pxpermm;
        newtrack.y_mm          = newtrack.y/newtrack.pxpermm;
        newtrack.a_mm          = newtrack.a/newtrack.pxpermm;
        newtrack.b_mm          = newtrack.b/newtrack.pxpermm;
        newtrack.theta_mm      = newtrack.theta;

        trx(i) = newtrack;
    end

    % Save trx.mat.
    save(fullfile(expdir,'trx.mat'),'timestamps','trx')

    % Write perframe folder.
    perframedir = fullfile(expdir,'perframe');
    if ~exist(perframedir,'dir')
      [success1,msg1] = mkdir(perframedir);
      if ~success1
        disp(msg1);
        return;
      end
    end

    n_feat = size(feat.data,3);
    for i=1:n_feat
        data = cell(1,n_flies);
        for s=1:n_flies
            data{s} = feat.data(s,trx(s).firstframe:trx(s).endframe,i);
        end
        units2.numerator = cell(1,0);
        units2.denominator = cell(1,0);
        save(fullfile(perframedir,[feat.names{i}]),'data','units2')
    end

    % Copy/soft-link movie.
    dosoftlink = 1;
    inmoviefile = moviename;
    if dosoftlink
        if exist(moviefile,'file')
            delete(moviefile);
        end
        if ispc
            if exist([moviefile,'.lnk'],'file')
                delete([moviefile,'.lnk']);
            end
            cmd = sprintf('mkshortcut.vbs /target:"%s" /shortcut:"%s"',inmoviefile,moviefile);
            fprintf('\nMaking a Windows shortcut file at "%s" with target "%s"\n',inmoviefile,moviefile);
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
        disp(' ')
        disp('WARNING: Could not write softlink');
        %copyfile(inmoviefile,moviefile);
    end

catch
     disp('WARNING: Could not write JAABA folders');
end