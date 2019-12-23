%% Header
%
% Publish all matlab code in the repository
%
% V.Pettas/F.Costa
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019

%%
clc,clear all %#ok<*CLALL>

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath (genpath('..\Functions'))
addpath (genpath('..\Main'))
addpath (genpath('..\HelpfulStandAlone'))

files  = dir('..\Main\*.m');
files2 = dir('..\Functions\*.m');
files3 = dir('..\HelpfulStandAlone\*.m');

for i = 1:size(files,1)
    names_main = cellstr(files(i).name );
    publish([ names_main{1}],'evalCode',false);
end
for i = 1:size(files2,1)
    names_functions = cellstr(files2(i).name ) ;
    publish([names_functions{1}],'evalCode',false);
    
end
for i = 1:size(files3,1)
    names_help = cellstr(files3(i).name)  ;
    publish([ names_help{1}],'evalCode',false);
end





