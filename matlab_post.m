%% plates-shells visualization

clc
clear all
close all

%% main

nts = 20;
num_len = 50;
num_wid = 40;

F(nts) = struct('cdata',[],'colormap',[]);

filename = strcat('result', num2str(0), '.txt');
data = importdata(filename);

x = reshape(data(:,1), [num_len num_wid]);
x = transpose(x(:,1));
y = reshape(data(:,2), [num_len num_wid]);
y = y(1,:);
[X, Y]=meshgrid(x,y);

z = reshape(data(:,3), [num_len num_wid]);
Z = transpose(z);
mesh(X,Y,Z);

xlim([0 4]);
ylim([0 3]);
zlim([0 2]);

for i = 1:nts
    filename = strcat('result', num2str(i), '.txt');
    data = importdata(filename);
    
    x = reshape(data(:,1), [num_len num_wid]);
    x = transpose(x(:,1));
    y = reshape(data(:,2), [num_len num_wid]);
    y = y(1,:);
    [X, Y]=meshgrid(x,y);
    
    z = reshape(data(:,3), [num_len num_wid]);
    Z = transpose(z);
    mesh(X,Y,Z);    
    
    xlim([0 4]);
    ylim([0 3]);
    zlim([0 2]);
    drawnow;
    F(i) = getframe;
    
    
end