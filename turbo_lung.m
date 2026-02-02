function turbo_lung(folder, subject,inpfunc,vox)
% TURBO_LUNG  Perform lung kinetic modeling from PET/CT data for
%             integration in Turku's TURBO pipeline
%             [15O]H2O with PA or RV input.
%
%   turbo_lung(folder, subject)
%
%   This function performs:
%     - Unzipping NIfTI files
%     - Loading PET and CT segmentation volumes
%     - Creating anatomical masks (PA, RV, lung)
%     - Extracting TACs from regions
%     - Generating PMOD-compatible TAC files
%     - ROI-level kinetic modeling (fit_h2o)
%     - Voxelwise modeling (turbo_vox_modelling)
%     - Writing parametric images and summary resultswait<
%
%   INPUTS:
%       folder  : base directory for subject data
%       subject : subject ID (string)
%       inpfunc : input function, 0 for PA, 1 for RV
%
%   REQUIRES:
%       fitdelay.m, fit_h2o.m, turbo_vox_modelling.m, get_germano.m
%
%   Author: Lars Tolbod
%   Date:   Nov 2025
% -------------------------------------------------------------------------

close all
%% === Unzip NIfTI files ===
disp('Unzipping NIfTI files...')
A = dir(fullfile([folder '/nifti/' subject], '**', '*.nii.gz'));
for i = 1:length(A)
    gunzip(fullfile(A(i).folder, A(i).name));
end

%% === Load anatomical and PET data ===
disp('Loading data...')
segHeartV = medicalVolume([folder '/nifti/' subject '/anat/' subject '_ct_segmentations_heart_coregistered.nii.gz']);
segV      = medicalVolume([folder '/nifti/' subject '/anat/' subject '_ct_segmentations_coregistered.nii.gz']);
petV      = medicalVolume([folder '/nifti/' subject '/pet/' subject '_pet.nii.gz']);

% Load options
turbo_options = load([folder '/nifti/' subject '/options/' subject '_turbo_options.mat']);
model_options = load([folder '/nifti/' subject '/options/' subject '_modelling_options.mat']);
study_options = load([folder '/nifti/' subject '/options/' subject '_study_specs.mat']);

% Identify lung ROI labels
lung_label = find(contains(turbo_options.roinames, 'lung'));

%% === Create masks ===
disp('Creating masks...')

petVsum = petV;
petVsum.Voxels = sum(petV.Voxels, 4);

se = strel('cube', 3);

% Pulmonary artery (PA)
maskPA = segHeartV.Voxels == 7;
maskPA = logical(imerode(maskPA, se));
maskPA = logical(imerode(maskPA, se));

% Right ventricle (RV)
maskRV = segHeartV.Voxels == 5;
thresRV=prctile(petVsum.Voxels(maskRV),75,"all");
tmp=petVsum.Voxels.*maskRV;
maskRV=tmp>thresRV;
maskRV = imclose(maskRV, se);
maskRV = imfill(maskRV, 'holes');
maskRV = logical(imerode(maskRV, se));

% Lung masks
maskTotal = false(size(segV.Voxels));
maskLungIdx_ = [];
for i = 1:length(lung_label)
    maskLung{i} = segV.Voxels == lung_label(i);
    maskTotal = maskTotal | maskLung{i};
    maskLung{i} = logical(imerode(maskLung{i}, se));
    maskLungIdx{i} = find(maskLung{i});    
    maskLungIdx_ = [maskLungIdx_; maskLungIdx{i}];
end

% Close and fill total lung mask
se = strel('sphere', 3);
maskTotal = imclose(maskTotal, se);
maskTotal = imfill(maskTotal, 'holes');
maskTotalidx=find(maskTotal);

%% === Compute regional TACs ===
disp('Extracting TACs...')

nFrames = size(petV.Voxels, 4);
for i = 1:nFrames
    fr = petV.Voxels(:,:,:,i);
    tacPA(i)        = mean(fr(maskPA), 'all');
    tacRV(i)        = mean(fr(maskRV), 'all');
    tacTotalLung(i) = mean(fr(maskTotal), 'all');
    for j = 1:length(lung_label)
        tacLung(i,j) = mean(fr(maskLung{j}), 'all');
    end
end

%% === Plot TACs ===
try 
    T = readtable([folder '/nifti/' subject '/input/' subject '_IDIFaorta_mc_rigid.txt']);
catch
    T = readtable([folder '/nifti/' subject '/input/' subject '_IDIFaorta_no_mc.txt']);
end
figure;
plot(T.Var1, tacPA, 'r',T.Var1, tacRV, 'g', T.Var1, tacTotalLung, 'b');
hold on
for i = 1:length(lung_label)
    plot(T.Var1, tacLung(:,i))
end
xlim([0 1.5])
roinames = strrep(turbo_options.roinames(lung_label), '_', ' ');
legend(["PA","RV", "Total Lung", roinames'])


if inpfunc==0
    tacIF=tacPA;
    nameIF='IDIFPA';
    disp('Using PA input')
elseif inpfunc==1
    tacIF=tacRV;
    nameIF='IDIFRV';
    disp('Using RV input')
end
if ~isfolder([folder '/nifti/' subject '/results/results_no_mc/roi-modelling/fit_h2o'])
    mkdir([folder '/nifti/' subject '/results/results_no_mc/roi-modelling/fit_h2o'])
end
saveas(gcf, [folder '/nifti/' subject '/results/results_no_mc/roi-modelling/fit_h2o/' subject '_' nameIF '_lung_tacs.png'])
%% === Generate PMOD-compatible TAC files ===
disp('Saving PMOD TAC files...')
time_min   = T.Var1;
wb_kBq     = tacIF / 1000;
pl_kBq     = tacIF / 1000;
time_start = study_options.frames(:,1);
time_end   = study_options.frames(:,2);
tissue1    = tacTotalLung / 1000;
tissue2    = tacLung(:,1) / 1000;

% Whole Blood file
blood_data = table(time_min, wb_kBq', 'VariableNames', {'time[min]', 'plasma[kBq/cc]'});
blood_filename = [folder '/nifti/' subject '/input/' subject '_pmod_RV_WholeBlood_TAC.txt'];
fid = fopen(blood_filename, 'w'); fprintf(fid, '#PMOD_BLD/n'); fclose(fid);
writetable(blood_data, blood_filename, 'Delimiter', 'tab');
fprintf('✅ Saved: %s/n', blood_filename);

% Tissue TAC file
tissue_data = table(time_start, time_end, tissue1', tissue2, ...
    'VariableNames', {'start[min]', 'end[min]', 'Tissue1[kBq/cc]', 'Tissue2[kBq/cc]'});
tissue_filename = [folder '/nifti/' subject '/input/' subject '_pmod_Lung_Tissue_TAC.txt'];
fid = fopen(tissue_filename, 'w'); fprintf(fid, '#PMOD_TAC/n'); fclose(fid);
writetable(tissue_data, tissue_filename, 'Delimiter', 'tab');
fprintf('✅ Saved: %s/n', tissue_filename);

%% === ROI-level modeling ===
disp('Performing ROI-level modeling...')
inputdata = double([mean(study_options.frames,2)*60 tacIF']);
frametimes = study_options.frames*60;
tac = double([tacLung tacTotalLung']);
delayfitmodel = 1;
verbose = 5;

figure; tiledlayout('flow')
for i = 1:size(tac,2)
    % Delay fitting
    [~, dly_corrected_inputdata, fitdelay_result(i), ~, ~, ~] = fitdelay( ...
        frametimes, tac(:,i)', inputdata, delayfitmodel, 'verbose', verbose);

    % Model fitting
    f = 0; fa = 1;
    n_iter = model_options.n_iter;
    lb = model_options.lb(1,:);
    ub = model_options.ub(1,:);
    [modelfit, x_optim, flow, resnorm, residual] = fit_h2o( ...
        frametimes, tac(:,i), ...
        dly_corrected_inputdata(:,1), dly_corrected_inputdata(:,2), ...
        dly_corrected_inputdata(:,2), f, fa, n_iter, lb, ub);

    % Extract parameters
    K1(i)     = x_optim(1);
    k2(i)     = x_optim(1)/x_optim(2);
    pWater(i) = x_optim(2);
    Va(i)     = x_optim(3);
    delay(i)  = x_optim(4);
    flow(i)   = (1 - Va(i))*K1(i);

    % Scale to standard units
    K1(i)   = K1(i)*100*60;
    k2(i)   = k2(i)*60;
    flow(i) = flow(i)*100*60;
    Va(i)   = Va(i)*100;

    % Plot model fit
    nexttile
    plot(inputdata(:,1), modelfit, 'b', ...
         inputdata(:,1), tac(:,i), 'x', ...
         inputdata(:,1), inputdata(:,2)*Va(i)/100, 'r');
    if i <= length(lung_label)
        title({roinames(i), ...
            sprintf('K1=%.3g  k2=%.3g  Va=%.2f  Delay=%.2f', ...
            K1(i), k2(i), Va(i), fitdelay_result(i)+delay(i))})
    else
        title({'Total Lung', ...
            sprintf('K1=%.3g  k2=%.3g  Va=%.2f  Delay=%.2f', ...
            K1(i), k2(i), Va(i), fitdelay_result(i)+delay(i))})
    end

    % Goodness of fit
    r = corrcoef(tac(:,i), modelfit);
    R_2(i) = r(1,2)^2;
end

saveas(gcf, [folder '/nifti/' subject '/results/results_no_mc/roi-modelling/fit_h2o/' subject '_' nameIF '_lung_ROIfit.png'])


%% === Voxelwise modeling ===
vox_results_dir = [folder '/nifti/' subject '/results/results_no_mc/vox-modelling/imgh2obf/'];
if ~isfile([vox_results_dir subject '_k1_imgh2obf_' nameIF '_lung_from_smoothed.nii']) || vox
disp('Performing voxelwise modeling...')
vox_petfile = [folder '/nifti/' subject '/pet/' subject '_pet_smoothed.nii'];
if ~isfile(vox_petfile)
    petfile=[folder '/nifti/' subject '/pet/' subject '_pet.nii'];
    smoothPET(petfile,subject)
end
frametimes  = study_options.frames*60;
inputdata   = double([mean(study_options.frames,2)*60 tacIF']);
maskMat     = load([folder '/nifti/' subject '/pet/' subject '_pet_mask.mat']);
fit_len_timeunit = 'sec';
vox_par_est = [];

%for i = 1:length(maskLungIdx)
    [vox_par_est, vox_par_est_names, vox_par_est_units] = turbo_vox_modelling( ...
        model_options.VOXmodel, vox_petfile, frametimes, inputdata, inputdata, ...
        mask_idx = maskTotalidx, ...
        verbose = verbose, ...
        fitdelay_start = fitdelay_result(i)+delay(i), ...
        fitdelay_end = fitdelay_result(i)+delay(i), ...
        timestep = 1, ...
        fit_va = model_options.fit_va, ...
        calculate_R2_img = model_options.calculate_R2_img, ...
        fit_length = model_options.vox_fit_length, ...
        fit_length_timeunit = fit_len_timeunit, ...
        start_time = model_options.vox_start_time, ...
        end_time = model_options.vox_end_time, ...
        use_parfor = turbo_options.use_parfor, ...
        weight = study_options.weight, ...
        dose = study_options.dose, ...
        LC = model_options.vox_LC);
    %vox_par_est = [vox_par_est; vox_par_est_];
%end

%% === Write parametric images ===
disp('Writing parametric images...')
% Gaussian smoothing for MIP visualization
voxelSpacing = petV.VoxelSpacing;
sigma_mm = 5;
sigma_voxels = sigma_mm ./ voxelSpacing;
if ~isfolder(vox_results_dir); mkdir(vox_results_dir); end
pethdr = spm_vol(vox_petfile); pethdr = pethdr(1);
vox_img_names = cell(length(vox_par_est_names), 1);

for i = 1:length(vox_par_est_names)
    if ~(isnan(model_options.filter_size) || model_options.filter_size == 0)
        pethdr.fname   = [vox_results_dir subject '_' lower(char(vox_par_est_names(i))) '_' model_options.VOXmodel '_' nameIF '_lung_from_smoothed.nii'];
    else
        pethdr.fname   = [vox_results_dir subject '_' model_options.VOXmodel '_' lower(char(vox_par_est_names(i))) '.nii'];
    end
    if strcmp(model_options.VOXmodel, 'imgh2o')
        if model_options.fit_va > 0
            pethdr.fname = strrep(pethdr.fname, '.nii', '_fitted_va.nii');
        else
            pethdr.fname = strrep(pethdr.fname, '.nii', '_fixed_va0.nii');
        end
    end
    pethdr.descrip = [subject ' ' char(vox_par_est_names(i)) '-parametric image'];
    vox_img_names{i} = pethdr.fname;
    temp_param_img = zeros(pethdr.dim);
    temp_param_img(maskTotalidx) = vox_par_est(:,i);
    temp_param_img=imgaussfilt3(temp_param_img, sigma_voxels);
    spm_write_vol(pethdr, temp_param_img);
end
end
%% === Visualization & summary ===
disp('Creating MIP visualizations...')
Vva    = medicalVolume([vox_results_dir subject '_va_imgh2obf_' nameIF '_lung_from_smoothed.nii']);
VK1    = medicalVolume([vox_results_dir subject '_k1_imgh2obf_' nameIF '_lung_from_smoothed.nii']);
Vk2    = medicalVolume([vox_results_dir subject '_k2_imgh2obf_' nameIF '_lung_from_smoothed.nii']);
Vdelay = medicalVolume([vox_results_dir subject '_delay_imgh2obf_' nameIF '_lung_from_smoothed.nii']);
Vflow  = medicalVolume([vox_results_dir subject '_flow_imgh2obf_' nameIF '_lung_from_smoothed.nii']);

% % Gaussian smoothing for MIP visualization
% voxelSpacing = VK1.VoxelSpacing;
% sigma_mm = 5;
% sigma_voxels = sigma_mm ./ voxelSpacing;
% 
% filteredVk1    = imgaussfilt3(VK1.Voxels, sigma_voxels);
% filteredVk2    = imgaussfilt3(Vk2.Voxels, sigma_voxels);
% filteredVva    = imgaussfilt3(Vva.Voxels, sigma_voxels);
% filteredVdelay = imgaussfilt3(Vdelay.Voxels, sigma_voxels);
% filteredVflow  = imgaussfilt3(Vflow.Voxels, sigma_voxels);
VK1vox=VK1.Voxels; Vk2vox=Vk2.Voxels; Vvavox=Vva.Voxels; 
Vdelayvox=Vdelay.Voxels; Vflowvox=Vflow.Voxels;
% Display MIPs
fh = figure; fh.WindowState = 'fullscreen'; tiledlayout(2,2)
titles = {'K1','k2','Va','Delay'};
vars   = {VK1vox, Vk2vox, Vvavox, Vdelayvox};
ranges = {[0 6],[0 20],[0 100],[0 5]};

for i = 1:4
    nexttile
    mip_coronal = rot90(squeeze(max(vars{i}, [], 2)));
    if i==1
        mip_coronal=mip_coronal./100;
    end
    imshow(mip_coronal, ranges{i}, 'colormap', get_germano)
    title(titles{i}); colorbar
end
saveas(gcf, [vox_results_dir subject '_' nameIF 'fixed_delay_lung_MiP.png'])

%% === Save ROI summary ===
disp('Saving ROI summary...')
for i = 1:length(maskLungIdx)
    results(i,:) = [ ...
        mean(VK1vox(maskLungIdx{i})), ...
        mean(Vk2vox(maskLungIdx{i})), ...
        mean(Vvavox(maskLungIdx{i})), ...
        mean(Vflowvox(maskLungIdx{i}))];
end

results(1+length(maskLungIdx),:) = [ ...
        mean(VK1vox(maskTotalidx)), ...
        mean(Vk2vox(maskTotalidx)), ...
        mean(Vvavox(maskTotalidx)), ...
        mean(Vflowvox(maskTotalidx))];


delay_tmp  = fitdelay_result + delay;
resultsROI = [K1(1:end)' k2(1:end)' Va(1:end)' R_2(1:end)' delay_tmp(1:end)'];
results    = [[roinames; 'Total Lung'] results resultsROI];

T = array2table(results, 'VariableNames', ...
    {'Segment','K1 (param)','k2 (param)','Va (param)','flow (param)', ...
     'K1 (ROI)','k2 (ROI)','Va (ROI)','R2 (ROI)','delay (ROI)'});
writetable(T, [vox_results_dir subject '_' nameIF '_results.xlsx']);
disp('✅ All processing completed.')


%% === Zip NIfTI files ===
disp('Zipping NIfTI files...')
A = dir(fullfile([folder '/nifti/' subject], '**', '*.nii'));
for i = 1:length(A)
    gzip(fullfile(A(i).folder, A(i).name));
    delete(fullfile(A(i).folder, A(i).name));
end



end


function smoothPET(petfile, subject)
        [pet_folder,pet_filename,~]=fileparts(petfile);
        pet_img = nifti(petfile);
        modelling_options.filter_size=3; verbose=5;
        pethdr = spm_vol(petfile);
        framenr = length(pethdr);
        pethdr = pethdr(1);
    
        pixdim  = abs([pethdr.mat(1,1)  pethdr.mat(2,2)  pethdr.mat(3,3)]);
        % dimfilt = round(11 ./ pixdim);
        fwhm    = modelling_options.filter_size ./ pixdim;
    
        % fwhm    = 3 ./ pixdim;
%         fwhm    = 5 ./ pixdim;
        sd      = fwhm/sqrt(8*log(2));
    
        temp_frame_paths = cell(framenr, 1);
    
        % Print progress
        if verbose > 0
            prog = 0.1;
            prog_print = sprintf("Smoothing dynamic image...");
            prog_print = pad(prog_print, 50, 'right');
    
            fprintf(prog_print); 
        end
        smooth_folder = [pet_folder filesep 'temp_smooth' filesep];
        if exist(smooth_folder, 'dir')
            delete([smooth_folder '*'])
        else
            mkdir(smooth_folder)
        end
        % Smooth frame by frame.
        
        for f = 1:framenr

            % smoothed_frame_data = PSVsmooth_3d(pet_img.dat(:,:,:,f), sd, dimfilt);
            smoothed_frame_data = imgaussfilt3(pet_img.dat(:,:,:,f), sd);

            % Save temporary smoothed frame.
            pethdr = spm_vol(petfile);
            pethdr = pethdr(1);
            temp_frame_file = [smooth_folder subject '_temp_smoothed_frame' int2str(f) '.nii'];
            pethdr.fname = temp_frame_file;
            temp_frame_paths{f} = temp_frame_file;
            spm_write_vol(pethdr, smoothed_frame_data);

        end
        % Clear memory

    
        petfile_smoothed = [pet_folder filesep pet_filename '_smoothed.nii'];
    
        spm_file_merge(temp_frame_paths, petfile_smoothed);
        delete([smooth_folder '*'])
end