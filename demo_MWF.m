addpath('.~/coreFunction/')
addpath('.~/MEDI_toolbox/')
addpath('.~/STISuite_V3.0/')
%
dicomFolderName = '~/mGREdata/';
[complex,voxelSize,matrix_size,CF,delta_TE,TE,B0_dir,B0] = Read_DICOM_HW(dicomFolderName);
%[mag,phase] = bipolarCorr(complex,matrix_size,voxelSize,TE);
[mag,phase] = bipolarCorrKspace(complex,matrix_size,voxelSize,TE);
mag  = double(mag);
phase = double(phase);
sz  =size(mag);
%% NLM based denoising
deMag = NLMDenoise(mag);

%% mask generated
fraction = 0.1;
Mask = maskGen(deMag,voxelSize,fraction);

%% Voxel spread Function to correct B0
MyInfo.FirstTE = TE(1);
MyInfo.Vox = voxelSize.*1e-3;
MyInfo.EchoSpacing = delta_TE;
MyInfo.EchoIndex = [1,3];
[magB0Corr,phaseB0Corr] = VSFB0Correction(deMag,phase,MyInfo);

%% if ram is small, using LFG
%%LFG correct B0
%MyInfo.Mask = Mask(:,:,:);
%MyInfo.Method = 1;	% Du method;MyInfo.Method = 2;	% other method
%MyInfo.FirstTE = 2.41e-3;
%MyInfo.EchoSpacing = 0.0014;
%MyInfo.Vox = [1.0 1.0 1.0] * 1e-3;
%MyInfo.BipolarFlag = true;
%[magB0Corr, Gp, Gv, Gs] = LFGB0Correction(mag,phase,MyInfo);

%% STISuIte based phase synthesising

padsize=[12 12 12];
smvsize = 12;
TempPhase = phaseSynth(phase,Mask,voxelSize,B0_dir,TE,B0,smvSize,padSize);

complexData = double(magB0Corr).*exp(-1i*double(TempPhase));%复数形式为乘以指数
%% BG initial estimating
BGInit = backGroundInit(complexData);

%% Fitting MWF using NLLS
% based on field strength and samples
%3T invivo
upperBound   = [2.0 0.024 2.0 0.150 2.0 0.150 75 25 25 0];
lowerBound   = [0.0 0.003 0.0 0.024 0.0 0.024 -75 -25 -25 0];
initialValue = [0.1 0.010 0.6 0.064 0.3 0.048 0 0 0 0];
MWF = myelinWaterFitting(complexData,Mask,BGInit,upperBound,lowerBound,initialValue);
%% save MWF.nii 
nii = make_nii(imrotate(MWF,90),voxelSize,[],64);
save_nii(nii,'MWF.nii')