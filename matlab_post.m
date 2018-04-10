%% plates-shells visualization

clc
clear all
close all

%% main

nts = 6;

for i = 1:nts
    filename = strcat('result', num2str(i), '.txt');
    data = importdata(filename);
    
    x = data(:,1);
    y = data(:,2);
    z = data(:,3);
    
    scatter3(x,y,z);
end