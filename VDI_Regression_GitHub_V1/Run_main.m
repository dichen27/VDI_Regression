clear;
clc;
close;

% Version 1.0
% Email:dichen727@gmail.com

%% add the  VDI functions
addpath(genpath('C:/VDI_Regression_GitHub_V1/cd_function/'));
working_dir='C:/VDI_Regression_GitHub_V1/';
cd(working_dir);

%% load example data
load('./data_sample.mat')

%% VDI for MRI type-1

VDI_MRI_1=cd_VDI_score(MRI_1); 

%% VDI for MRI type-2

VDI_MRI_2=cd_VDI_score(MRI_2);

%% Combine the MRI-1 and MRI-2

MRI_mix=[MRI_1,MRI_2];% Combine MRI-1 and MRI-2
VDI_mix=[VDI_MRI_1;VDI_MRI_2];% Combine VDI-1 and VDI-2

%% VDI regression for phenotype_1

brain_id=[ID_intersect,MRI_mix];% ID and combined MRI
score_id=[ID_intersect,phenotype_1];% ID and phenotype_1
t_1=cd_VDI_t_new(brain_id,score_id,cov);% for whole-brain t-map
t_square=t_1.*t_1;% t^2
sample_size_1=length(ID_intersect);% sample_size_1 is sample size for phenotype_1

% VDI regression of phenotype_1
% h_1: the BAV for phenotype_1
% intercept_1: the intercept for phenotype_1
[h_1,intercept_1,h_foundation,intercept_foundation]=cd_VDI_regression(VDI_mix,t_square,sample_size_1,1)

%% VDI regression for phenotype_2
brain_id=[ID_intersect,MRI_mix];% ID and combined MRI
score_id=[ID_intersect,phenotype_2]; % ID and phenotype_2
t_2=cd_VDI_t_new(brain_id,score_id,cov); % for whole-brain t-map
t_square=t_2.*t_2;% t^2
sample_size_2=length(ID_intersect);% sample_size_2 is sample size for phenotype_2


% VDI regression of phenotype_2
% h_2: the BAV for phenotype_2
% intercept_2: the intercept for phenotype_2
[h_2,intercept_2,h_2_foundation,intercept_2_foundation]=cd_VDI_regression(VDI_mix,t_square,sample_size_2,1)

%% interaction term

t_12=t_1.*t_2;
phenotype_r=corr(phenotype_1,phenotype_2);% behavioral correlation between phenotype_1 and phenotype_2
sample_size_intersect=length(ID_intersect);% sample_size_intersect is intersect sample size between phenotype_1 and phenotype_2

%% Neuroimaging correlation

% caculate the Covariance
[Covariance,intercept]=cd_VDI_corr(VDI_mix,t_12,sample_size_1,sample_size_2,h_foundation,h_2_foundation,1,1,sample_size_intersect,phenotype_r);

% Neuroimaging correlation r
r_brain=Covariance/sqrt(h_1*h_2)

%% finished







