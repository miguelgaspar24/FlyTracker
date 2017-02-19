
% Chops up excess data columns from feat-mat files that might have grown
% to 17 columns but ccomputed the features incorrectly and need to be
% re-run.

filepath = 'C:\Users\Miguel\Desktop\Cecilia\video128_2017-01-12T15_44_56\video128_2017-01-12T15_44_56-feat.mat';

names = feat.names(1:9);
units = feat.units(1:9);
data = feat.data(:,:,1:9);

feat.names = names
feat.units = units
feat.data = data

save(filepath, 'feat')

clearvars('feat')
disp('The "feat" variable has been cleared.')