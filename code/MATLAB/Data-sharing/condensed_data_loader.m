%% Header
% Condensed_Data_Loader.m
% Script/functions designed to load a single experiments data for external
% analysis

clear

% UI pop-up to find .mat file
[data_filename, data_path] = uigetfile('*.mat','Select data file (must be .mat)');
load([data_path data_filename])

clear data_filename data_path