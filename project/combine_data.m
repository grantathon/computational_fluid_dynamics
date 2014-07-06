function [] = combine_data(folder_name, variable, mean, stdDev)
close all

% obtain the name of all *.mc files in the folder
listing = dir(['/home/ayman/CSE/CFD/computational_fluid_dynamics/project/' folder_name '/*.mc']);

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

figure; hist(rv,20)
title([variable ': ' num2str(mean) '+/-' num2str(stdDev) ' N: ' num2str(files)])
xlabel('Reynolds #')

figure; hist(x,20)
title([variable ': ' num2str(mean) '+/-' num2str(stdDev) ' N: ' num2str(files)])
xlabel('Separation Point x')

figure; hist(t,20)
title([variable ': ' num2str(mean) '+/-' num2str(stdDev) ' N: ' num2str(files)])
xlabel('Steady state time t')