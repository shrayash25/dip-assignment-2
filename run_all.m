% EEE F435 Digital Image Processing - Assignment 2
% Run from MATLAB:  run_all
% Run from a shell:  matlab -batch "run_all"
% All figures go to results/, console output goes to results/log.txt

clear; close all; clc;
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));

ID = '2023A3PS0361P';
A = str2double(ID(4));    % 3
B = str2double(ID(10));   % 3

out = fullfile(root, 'results');
if ~exist(out, 'dir'), mkdir(out); end
if exist(fullfile(out, 'log.txt'), 'file'), delete(fullfile(out, 'log.txt')); end
diary(fullfile(out, 'log.txt'));
set(0, 'DefaultFigureVisible', 'off');
rng(0);

fprintf('BITS ID %s -> A = %d, B = %d\n', ID, A, B);

I1 = load_gray(fullfile(root, 'images', 'image1.jpg'), 512);
I2 = load_gray(fullfile(root, 'images', 'image2.jpg'), 512);
imgs = {I1, I2};
tags = {'img1', 'img2'};

for k = 1:2
    fprintf('\n===== %s =====\n', tags{k});
    task1_properties(imgs{k}, tags{k}, A, B, out);
    task2_filters(imgs{k}, tags{k}, out);
    task2_bands(imgs{k}, tags{k}, out);
    task2_periodic_noise(imgs{k}, tags{k}, A, B, out);
    task2_compression(imgs{k}, tags{k}, out);
    extras(imgs{k}, tags{k}, out);
end

task2_mag_phase(I1, I2, out);
task2_synthetic(out);

set(0, 'DefaultFigureVisible', 'on');
fprintf('\nDone. Figures saved in %s\n', out);
diary off;
