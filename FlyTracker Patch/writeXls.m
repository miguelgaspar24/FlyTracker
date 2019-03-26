
function writeXls(savename,data,names,trk,n_flies,flies,chamber,s,pair,counts)    
    warning off
    [~,msg] = xlswrite([savename '.xlsx'],ones(10));
    use_xls = ~(numel(msg.message)>8 && ...
                strcmp(msg.message(1:9),'Could not'));
    
    n_frames = size(data,2);
    n_feats = size(data,3);
    
    if isfield(trk,'flies_in_chamber')
        if use_xls
            delete([savename '.xlsx'])
            % Write to xls.
            if n_flies <= 2
                for i=1:n_flies
                    sheet = [chamber(2:end), '_subject', num2str(flies(i))];
                    dat = data(i,:,:);
                    dat = reshape(dat,n_frames,n_feats);
                    xlswrite([savename, chamber, '.xlsx'],names,sheet);
                    xlswrite([savename, chamber, '.xlsx'],dat,sheet,'A2');
                end
            elseif n_flies > 2
                sheet_no = 1;
                for i=[s, pair]
                    sheet = [chamber(2:end), '_subject', num2str(i)];
                    dat = data(sheet_no,:,:);
                    dat = reshape(dat,n_frames,n_feats);
                    xlswrite([savename, chamber, counts, '.xlsx'],names,sheet);
                    xlswrite([savename, chamber, counts, '.xlsx'],dat,sheet,'A2');
                    sheet_no = sheet_no + 1;
                end
            end   
        else    
            delete([savename '.csv'])
            % Write to csv folder.
            if ~exist(savename,'dir')
                mkdir(savename);
            end
            if n_flies <= 2
                for i=1:n_flies
                    dat = data(i,:,:);
                    dat = reshape(dat,n_frames,n_feats);
                    filename = fullfile(savename, [chamber(2:end), '_subject', num2str(flies(i)), '.csv']);
                    fid = fopen(filename,'w');
                    fprintf(fid, '%s,', names{1:end-1});
                    fprintf(fid, '%s\n', names{end});
                    fclose(fid);
                    dlmwrite(filename,dat,'-append');
                end
            elseif n_flies > 2
                sheet_no = 1;
                for i=[s, pair]
                    dat = data(sheet_no,:,:);
                    dat = reshape(dat,n_frames,n_feats);
                    filename = fullfile(savename, [chamber(2:end), '_pair', counts, '_subject', num2str(i), '.csv']);
                    fid = fopen(filename,'w');
                    fprintf(fid, '%s,', names{1:end-1});
                    fprintf(fid, '%s\n', names{end});
                    fclose(fid);
                    dlmwrite(filename,dat,'-append');
                    sheet_no = sheet_no + 1;
                end
            end
        end


    elseif ~isfield(trk,'flies_in_chamber')
        if use_xls
            delete([savename '.xls'])
            % Write to xls.
            if n_flies <= 2
                for i=1:n_flies
                    sheet = ['subject', num2str(i)];
                    dat = data(i,:,:);
                    dat = reshape(dat,n_frames,n_feats);
                    xlswrite([savename, '.xlsx'],names,sheet);
                    xlswrite([savename, '.xlsx'],dat,sheet,'A2');
                end
            elseif n_flies > 2
                sheet_no = 1;
                for i=[s, pair]
                    sheet = ['subject', num2str(i)];
                    dat = data(sheet_no,:,:);
                    dat = reshape(dat,n_frames,n_feats);
                    xlswrite([savename, counts, '.xlsx'],names,sheet);
                    xlswrite([savename, counts, '.xlsx'],dat,sheet,'A2');
                    sheet_no = sheet_no + 1;
                end
            end   
        else    
            delete([savename '.csv'])
            % Write to csv folder.
            if ~exist(savename,'dir')
                mkdir(savename);
            end
            if n_flies <= 2
                for i=1:n_flies
                    dat = data(i,:,:);
                    dat = reshape(dat,n_frames,n_feats);
                    filename = fullfile(savename, ['subject', num2str(i), '.csv']);
                    fid = fopen(filename,'w');
                    fprintf(fid, '%s,', names{1:end-1});
                    fprintf(fid, '%s\n', names{end});
                    fclose(fid);
                    dlmwrite(filename,dat,'-append');
                end
            elseif n_flies > 2
                sheet_no = 1;
                for i=[s, pair]
                    dat = data(sheet_no,:,:);
                    dat = reshape(dat,n_frames,n_feats);
                    filename = fullfile(savename, ['pair', counts, '_subject', num2str(i), '.csv']);
                    fid = fopen(filename,'w');
                    fprintf(fid, '%s,', names{1:end-1});
                    fprintf(fid, '%s\n', names{end});
                    fclose(fid);
                    dlmwrite(filename,dat,'-append');
                    sheet_no = sheet_no + 1;
                end
            end
        end
    end
    warning on
end
