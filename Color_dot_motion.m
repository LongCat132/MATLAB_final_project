%% add data folder
if ~exist('Color_data/') 
    mkdir Color_data
    mkdir Color_data/subinfoQ
    mkdir Color_data/pre_limits
end 
%% %%%%%%%% General Settings %%%%%%%%%%%%%%  
% initiation
clear; clc; close all;
addpath('__self_func__');
rng('shuffle'); % shuffle the randome number seed every time when matlab restarts
[SubID,SubName,SubAge,SubGender,Handedness] = greeting();
disp(['Hello, ' SubName '!']);
save(['Color_data/subinfo/' num2str(SubID) '_' SubName '.mat'],...
            'SubID','SubName','SubAge','SubGender','Handedness');        
%% Pre-block,testing the subject's speed threshold
ListenChar(2);
[JND, done] = Speed_Judge(SubID,SubName);
%% Main Procedure
% rng needs to be set in the main procedure?
Color_main(SubID,SubName,JND);
ListenChar(0);
