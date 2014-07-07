function [] = combine_data()
close all

% obtain the name of all *.mc files in the folder
folders = dir(pwd);
i_folder = [folders(:).isdir];
folder_names = {folders(i_folder).name}';
folder_names(ismember(folder_names,{'.','..'})) = [];

num_fold = numel(folder_names);

rv = [];
x = [];
t = [];
N = 0;
for (i=1:num_fold)
    mc_file_list = dir([pwd '/' folder_names{i,1} '/*.mc']);
    
    % collected data from one particular folder
    [rv_tmp, x_tmp, t_tmp, N_tmp] = collect_mc(mc_file_list, folder_names{i,1});
    
    % combine data from each folder
    rv = [rv rv_tmp];
    x = [x x_tmp];
    t = [t t_tmp];
    N = N + N_tmp;
    
    clear rv_tmp x_tmp t_tmp
end

figure;
subplot(1,3,1); hist(rv,20)
ylabel('samples'); xlabel('Reynolds #')

subplot(1,3,2); hist(x,20)
title(['viscosity, N: ' num2str(N)])
ylabel('samples'); xlabel('Separation Point x [m]')

subplot(1,3,3); hist(t,20)
ylabel('samples'); xlabel('Steady state time t [sec]')

end

function [rv, x, t, files] = collect_mc(listing, folder_name)
[files, ~] = size(listing);

rv = zeros(1,files);
x = zeros(1,files);
t = zeros(1,files);

% open each *.mc file and make x and t vectors
for (i = 1:files)
    file_id = fopen([folder_name '/' listing(i).name]);
    
    tmp = fscanf(file_id,'%f %f %f');
    rv(i) = tmp(1);
    x(i) = tmp(2);
    t(i) = tmp(3);
    
    fclose(file_id);
end

end