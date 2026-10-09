function userOptions = projectOptions_fsl_demoGG(subjectNames)
% projectOptions_fsl_demo is an options file for the FSL demo tutorial
%
% Roni Maimon 9-2017 
% Georgie 2023
%__________________________________________________________________________
% Copyright (C)

%% Project details
%Current folder
curFolder=fileparts(which('projectOptions_fsl_demo.m'));

% This name identifies a collection of files which all belong to the same run of a project.
userOptions.analysisName = 'DEMO_FSL'; % this is renamed in the code for demos 3-4.

% This is the root directory of the project.
userOptions.rootPath = [curFolder,filesep,'DEMO_FSL'];

userOptions.betaPath = [pwd,filesep,'demoTest',filesep,'[[subjectName]]',filesep,'[[betaIdentifier]]'];



%%%%%%%%%%%%%%%%%%%%%%%%
%% EXPERIMENTAL SETUP %%
%%%%%%%%%%%%%%%%%%%%%%%%

% The list of subjects to be included in the study.
userOptions.subjectNames = subjectNames;

%% FSL Specific parameters
userOptions.run_names = { ...
    '01','02','03','04', '05', '06','07','08','09','10'};
userOptions.featsPrefix = 'run-';
userOptions.featsSuffix = '';
% userOptions.featsPath = [curFolder,filesep,'ProjectData',filesep,'[[subjectName]]',...
%                         filesep,'func',filesep,'task-moon',filesep,'model',filesep,'RSA', ...
%                         filesep, '[[featPrefix]]','[[runName]]','[[featSuffix]]','.feat'];
userOptions.featsPath = [curFolder,filesep,'Pipeline1',filesep,'[[subjectName]]',...
                        filesep,'func',filesep,'task-moon',filesep,'model',filesep,'RSA', ...
                        filesep, '[[featPrefix]]','[[runName]]','[[featSuffix]]','.feat'];

userOptions.copes = num2cell(1:84);

%% End of FSL specific parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYSIS PREFERENCES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% First-order analysis
% Text lables which may be attached to the conditions for MDS plots.
userOptions.conditionLabels = cellstr(["n" + [1:24,91:96], "g" + [1:24] "c" + [1:24,91:96]]);
userOptions.useAlternativeConditionLabels = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FEATUERS OF INTEREST SELECTION OPTIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%% %% %% %% %%
	%% fMRI  %% Use these next three options if you're working in fMRI native space:
	%% %% %% %% %%
	
	% The path to a stereotypical mask data file is stored (not including subject-specific identifiers).
	% "[[subjectName]]" should be used as a placeholder to denote an entry in userOptions.subjectNames
	% "[[maskName]]" should be used as a placeholder to denote an entry in userOptions.maskNames
	%userOptions.maskPath = [curFolder,filesep,'ProjectData/derivatives',filesep,'[[subjectName]]',filesep,'masks',filesep,'[[maskName]].nii.gz'];
	userOptions.maskPath = [curFolder,filesep,'Pipeline1/derivatives',filesep,'[[subjectName]]',filesep,'masks',filesep,'[[maskName]].nii.gz'];
	
end%function
