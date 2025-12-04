function turbo_kidney_petmr(subject)
% TURBO_KIDNEY_PETMR  Complete PET/MR kidney workflow (IDIF + VOI + voxelwise modelling)
%
% DESCRIPTION
%   This function performs a full kidney PET/MR processing pipeline for a
%   given subject folder. Utilizes TURBO. The workflow includes:
%
%     1. Directory setup and DICOM → NIfTI conversion (PET & MR)
%     2. MR–PET resampling
%     3. Whole-body segmentation using TotalSegmentator (MR-based)
%     4. Aorta extraction and refinement → image-derived input function (IDIF)
%     5. Kidney mask extraction (L/R), registration, cortical shell estimation
%     6. TAC extraction for IDIF and kidney cortices
%     7. VOI-based kinetic modelling using the H2O model (K1, k2, Va, delay, R²)
%     8. Voxelwise parametric imaging using the turbo_vox_modelling framework
%     9. Saving NIfTI parameter maps and summary tables (.xlsx)
%
%   Eventhough TotalSegmentator is trained on T1w images, it works well on 
%   InPhase 3D LAVA Flex images (GE Signa PET/MR).  
%
% INPUT
%   subject  – string containing the subject ID, which must correspond to
%              subfolders in your data_path and dicom_path directories.
%
% OUTPUT
%   No direct output variables. The function writes:
%       - NIfTI files (resampled MR, masks, smoothed PET, parametric images)
%       - PNG figures (fits, kidney VOIs)
%       - Excel tables with VOI and voxelwise results
%       - PMOD-compatible blood input file (.bld)
%
% REQUIREMENTS / DEPENDENCIES
%   - MATLAB Image Processing Toolbox
%   - SPM12 (for NIfTI I/O)
%   - TotalSegmentator (installed and available in system PATH)
%   - Custom functions (must be available in MATLAB path):
%         nifti_conversion
%         refineAortaMask
%         findBoundingBox
%         greatestConnectedComponentOnly
%         fitdelay
%         fit_h2o
%         turbo_vox_modelling
%
% NOTES
%   - Expects global variables or workspace vars: data_path, dicom_path.
%   - Assumes PET data is decay-corrected.
%   - H2O model options defined via set_modelling_options() and set_turbo_options().
%   - Segmentation labels correspond to TotalSegmentator MR task output.
%
% AUTHOR
%   Lars Tolbod / AUH
%   Dec 2025
%
% -------------------------------------------------------------------------



subject_dir = sprintf('%s/%s',data_path,subject);
anat_dir = sprintf('%s/anat',subject_dir);
pet_dir = sprintf('%s/pet', subject_dir);
input_dir = sprintf('%s/input', subject_dir); 

if ~exist(subject_dir,'dir')
    mkdir(subject_dir)  % Main folder
end
if ~exist(anat_dir,'dir')
    mkdir(anat_dir)  % Main folder
end
if ~exist(pet_dir,'dir')
    mkdir(pet_dir)  % Main folder
end
if ~exist(input_dir,'dir')
    mkdir(input_dir)  % Main folder
end

pet_orig_file = [pet_dir filesep subject '_pet.nii'];
mr_orig_file  = [subject_dir filesep 'anat' filesep subject '_mr.nii'];
mr_orig_seg_file  = [subject_dir filesep 'anat' filesep 'resliced_' subject '_mr_segmentations.nii'];
PET_dicom_folder = [dicom_path filesep subject filesep 'PET'];
MR_dicom_folder  = [dicom_path filesep subject filesep 'MR'];

%try lowercase if uppercase folders don't exist
if ~exist(PET_dicom_folder,'dir')
    PET_dicom_folder = [dicom_path filesep subject filesep 'pet']; 
end 
if ~exist(MR_dicom_folder,'dir')
    MR_dicom_folder = [dicom_path filesep subject filesep 'mr']; 
end 


%Convert to NIFTI
if ~isfile(pet_orig_file)
    [~,n,~]=fileparts(pet_orig_file);
    nifti_conversion(PET_dicom_folder, pet_dir, n)
end
if ~isfile(mr_orig_file)
    [~,n,~]=fileparts(mr_orig_file);
    nifti_conversion(MR_dicom_folder, anat_dir, n)
end

%Resample MR
petV       = medicalVolume(pet_orig_file);
mrV = medicalVolume(mr_orig_file);
mrVcoreg = resample(mrV, petV.VolumeGeometry, ...
                            Method="linear", FillValue=0);
[p,n,e]=fileparts(mr_orig_file);
mr_resampled_file=[p filesep n '_resampled' e];
if ~isfile(mr_resampled_file)
    write(mrVcoreg, mr_resampled_file);
end

%Run totalsegmentator
mr_seg_nii = strrep(mr_resampled_file,'.nii','_segmentations.nii');
if ~isfile(mr_seg_nii)
    [~, seg_version] = system('TotalSegmentator --version');
    system(['TotalSegmentator -i ' mr_resampled_file ' -o ' mr_seg_nii ' --ml --task total_mr']);
end
%Read
segV       = medicalVolume(mr_seg_nii);
petVsum    = petV;
petVsum.Voxels=sum(petV.Voxels,4);
petVsum.Voxels(:,:,1:2)=zeros(petVsum.NumCoronalSlices,petVsum.NumSagittalSlices,2);
petVsum.Voxels(:,:,petVsum.NumTransverseSlices-1:petVsum.NumTransverseSlices)=zeros(petVsum.NumCoronalSlices,petVsum.NumSagittalSlices,2);
study_specs=get_studyspecs(pet_orig_file);
%% 
frametimes = study_specs.frames;
inputTimes = mean(study_specs.frames,2);
model_options = set_modelling_options(tracer='h2o',frametimes=frametimes);
turbo_options = set_turbo_options(tracer='h2o',frametimes=frametimes);
model_options.vox_fit_length=240;

%%
pixdim  = petV.VoxelSpacing;
fwhm    = model_options.filter_size ./ pixdim;
sd      = fwhm/sqrt(8*log(2));
pet_orig_smooth_nii = strrep(pet_orig_file,'.nii','_smooth.nii');
petVsmooth=petV;
if ~isfile(pet_orig_smooth_nii)
    petVsmooth.Voxels = single(zeros(size(petV.Voxels)));  % Initialize smoothed PET volume
    for f = 1:size(petV.Voxels, 4)
        petVsmooth.Voxels(:,:,:,f) = single(imgaussfilt3(petV.Voxels(:,:,:,f), sd));
        pethdr = spm_vol(pet_orig_file);
        pethdr = pethdr(1);
        temp_frame_file = strrep(pet_orig_file,'.nii',['_tmp' int2str(f) '_smooth.nii']);
        pethdr.fname = temp_frame_file;
        temp_frame_paths{f} = temp_frame_file;
        spm_write_vol(pethdr, petVsmooth.Voxels(:,:,:,f));
    end
    spm_file_merge(temp_frame_paths, pet_orig_smooth_nii);
else
    petVsmooth=medicalVolume(pet_orig_smooth_nii);
end
%% Extract IDIF

baseMask=false(size(segV.Voxels));
baseMask(segV.Voxels==23)=true;
seAorta = strel('cube',2);

aortaMask = refineAortaMask(baseMask, petVsum.Voxels, seAorta);
if showDebug
volshow(petVsum,'OverlayData',aortaMask,'OverlayRenderingStyle','LabelOverlay','RenderingStyle','SlicePlanes')
end
aortaMaskV=petVsum;
aortaMaskV.Voxels=int8(aortaMask);
maskFile=[subject_dir filesep 'input' filesep subject '_bloodmask.nii'];
write(aortaMaskV,maskFile)
%% ==== IDIF extraction (aorta) ====
numFrames=size(petV.Voxels,4);
disp('Extracting IDIFs...');
idif = zeros(numFrames, 1);
for f = 1:numFrames
    frame = petV.Voxels(:,:,:,f);
    idif(f) = mean(frame(aortaMask));
end
idifFile  = [subject_dir filesep 'input' filesep subject '_idif.bld'];
inputdataFull = double([inputTimes idif/1000]);
save_pmod_bld(inputdataFull, idifFile)
%% ==== Registration optimizer ====
[optimizer, metric] = imregconfig('multimodal');
optimizer.InitialRadius       = 0.001;
optimizer.GrowthFactor        = 1.005;
optimizer.MaximumIterations   = 400;
% ---- Automatic voxel spacing detection ----
% medicalVolume.VoxelSpacing is [dy dx dz] or similar; we only need
% the physical size in each dimension and sometimes the mean.
voxelSpacing = petV.VoxelSpacing(:)';  % [sY sX sZ] (mm)
voxelSizeMean = mean(voxelSpacing);    % mm, used for isotropic morphology
dz = voxelSpacing(3);                  % slice thickness in mm
%% ==== Process each kidney ====
disp('Extracting cortical masks and kidney TACs...');
kidney_label=[2 3]; numKidneys=length(kidney_label); 
roiNames={'spleen','kidney_right','kidney_left'};
cortexMaskFull=zeros(size(petVsum.Voxels));
for ki = 1:numKidneys

    roiName = roiNames{kidney_label(ki)};
    if ki==1
        isRightKidney = true;
    else
        isRightKidney = false;
    end

    % Kidney binary mask from segmentation
    voxelMask = (segV.Voxels == kidney_label(ki));

    % Bounding box in index space
    [x1,y1,z1, x2,y2,z2] = findBoundingBox(voxelMask);
    zBounds(:,ki) = [z1; z2];

    % Crop smoothed PET & kidney mask
    petCrop     = petVsmooth.Voxels(y1:y2, x1:x2, z1:z2, :);
    petMeanCrop = mean(petCrop, 4);
    maskCrop    = voxelMask(y1:y2, x1:x2, z1:z2);

    % Registration: kidney mask -> PET
    maskSmooth = imgaussfilt3(single(maskCrop), 2.0);
    tform = imregtform(maskSmooth, petMeanCrop, 'rigid', optimizer, metric);
    mask_reg = imwarp( ...
        maskCrop, tform, ...
        'OutputView', imref3d(size(petMeanCrop)), ...
        'InterpolationMethod','nearest');   % keep binary

    % ---- Morphology: cortex shell based on physical thickness ----
    cortexThickness_mm = 10;                        % physical thickness
    rVox = max(1, round(cortexThickness_mm / voxelSizeMean));
    seCortex = strel('sphere', rVox);

    M = imfill(mask_reg, 'holes');
    M = imclose(M, strel('sphere', 20));            % fill dents

    innerKidney  = imerode(M, seCortex);
    cortexShell  = M & ~innerKidney & mask_reg;

    % ---- PET-based thresholding for cortex ----
    maskedPET = petMeanCrop .* cortexShell;
    threshold = prctile(petMeanCrop(:), 75);        % high-uptake shell
    thresMask = maskedPET > threshold;

    finalMask = imerode(thresMask, strel('sphere', 1));

    % Store 4D PET crop from ORIGINAL PET volume (not smoothed)
    kidneyPET{ki}   = petV.Voxels(y1:y2, x1:x2, z1:z2, :);
    cortexMasks{ki} = finalMask;
    cortexMaskFull(y1:y2, x1:x2, z1:z2)=cortexMasks{ki}*ki;
    
end

if showDebug
   
    volshow(petVsum, 'OverlayData', cortexMaskFull, ...
        'OverlayRenderingStyle','LabelOverlay', ...
        'RenderingStyle','SlicePlanes','OverlayAlpha',0.7);
end
kidneyMaskV=petVsum;
kidneyMaskV.Voxels=int8(cortexMaskFull);
kidneymaskFile=[subject_dir filesep 'anat' filesep subject '_kidney_masks.nii'];
write(kidneyMaskV,kidneymaskFile)

%% ==== TAC extraction (kidneys) ====
disp('Extracting kidney TACs...');
tac = zeros(numFrames, numKidneys);

for ki = 1:numKidneys
    mask = cortexMasks{ki};
    maskIdx = mask(:);   % speed: linear indexing

    % reshape kidneyPET to 2D [nVox x nFrames] for vectorized mean
    vol4D = kidneyPET{ki};
    vol2D = reshape(vol4D, [], size(vol4D,4));   % [nVox x nFrames]
    tac(:,ki) = mean(vol2D(maskIdx,:), 1).';     % mean over voxels
end

figure;
tiledlayout(2,2);
nexttile
plot(inputTimes/60,idif/1000,'r',...
    inputTimes/60,tac(:,1)/1000,'b',...
    inputTimes/60,tac(:,2)/1000,'g')
xlabel('Time (min)')
ylabel('kBq/mL')
legend({'IDIF',strrep(roiNames{kidney_label(1)},'_',' '),strrep(roiNames{kidney_label(2)},'_',' ')})

idif = idif ./ 1000;
tac  = tac  ./ 1000;



%% Run VOI kinetic analysis


range=inputTimes<model_options.vox_fit_length;
tmp=tac(:,1);
tac_(:,1)=tmp(range);
tmp=tac(:,2);
tac_(:,2)=tmp(range);
tmp=frametimes(:,1);
frametimes_(:,1)=tmp(range);
tmp=frametimes(:,2);
frametimes_(:,2)=tmp(range);

nIF    = size(idif,2);
nParam = 5;      % [K1 k2 Va R2 delay]
resultsH = zeros(1, nIF*nParam);
resultsV = zeros(1, nIF*nParam);

for id = 1:nIF
    tacIF = idif(:,id);
    inputdata = double([inputTimes(range) tacIF(range)]);
    

    for ki = 1:numKidneys

        [~, corrected, fitdelay_result, ~, ~, ~] = fitdelay( ...
            frametimes_, tac_(:,ki)', inputdata, 1, 'verbose', 5);

        [modelfit, x_opt, ~, ~, ~] = fit_h2o( ...
            frametimes_, tac_(:,ki), ...
            corrected(:,1), corrected(:,2), corrected(:,2), ...
            0, 1, ...
            model_options.n_iter, ...
            model_options.lb(1,:), model_options.ub(1,:));

        % Extract parameters
        K1   = x_opt(1) * 100 * 60;
        k2   = (x_opt(1)/x_opt(2)) * 60;
        Va   = x_opt(3) * 100;
        dly  = fitdelay_result + x_opt(4);
        flow = (1 - x_opt(3)) * x_opt(1) * 100 * 60;

        r  = corrcoef(tac_(:,ki), modelfit);
        R2 = r(1,2)^2;

        % Plot model fit
        nexttile;
        plot(inputdata(:,1), modelfit, 'b', ...
             inputdata(:,1), tac_(:,ki), 'x', ...
             inputdata(:,1), inputdata(:,2)*Va/100, 'r');
        title({[strrep(roiNames{kidney_label(ki)},'_',' ')], ...
               sprintf('K_1=%.3g  k_2=%.3g  V_a=%.2f  Delay=%.2f  R^2=%.2f', ...
                       K1, k2, Va, dly, R2)});
        xlabel('Time (s)');
        ylabel('Activity');

        idx0 = (id-1)*nParam + 1;
        if ki == 1
            resultsH(1, idx0:idx0+4) = [K1 k2 Va R2 dly];
        else
            resultsV(1, idx0:idx0+4) = [K1 k2 Va R2 dly];
        end
    end
end

%% === Voxelwise modeling ===
disp('Performing voxelwise modeling...')
vox_petfile = pet_orig_smooth_nii;
maskMat     = [];
fit_len_timeunit = 'sec';
vox_par_est = [];
verbose = 5;
% Inputs
V = mrVcoreg.Voxels;            % 3-D volume
% 1) Pre-smooth to reduce noise
Vsm = imgaussfilt3(V, 1.0);    % sigma in voxels; tune 0.5-2
% Vsm is your 3-D volume (numeric)
mask = (Vsm > 50);  
% Optional: remove small islands and smooth mask
mask = bwareaopen(mask, 1000);    % remove components smaller than 1000 voxels (tune)
se = strel('sphere',2);           % structural element radius in voxels (tune)
mask = imclose(mask, se);         % close small holes/gaps
mask = imopen(mask, se);          % remove small protrusions
inputdata_=inputdataFull;
inputdata_(:,2)=inputdata_(:,2)*1000;
mask_idx=find(mask);
[vox_par_est_, vox_par_est_names, vox_par_est_units] = turbo_vox_modelling( ...
        model_options.VOXmodel, vox_petfile, frametimes, inputdata_, ...
        mask_idx = mask_idx, ...
        verbose = verbose, ...
        fitdelay_start = -5, ...
        fitdelay_end = 5, ...
        timestep = 1, ...
        fit_va = model_options.fit_va, ...
        calculate_R2_img = model_options.calculate_R2_img, ...
        fit_length = model_options.vox_fit_length, ...
        fit_length_timeunit = fit_len_timeunit, ...
        start_time = 0, ...
        end_time = 240/60, ...
        use_parfor = turbo_options.use_parfor, ...
        weight = 1, ...
        dose = study_specs.dose, ...
        LC = model_options.vox_LC);
vox_par_est = [vox_par_est; vox_par_est_];

%% === Write parametric images ===
disp('Writing parametric images...')
vox_results_dir = fullfile(subject_dir,'results','results_no_mc','vox-modelling','imgh2obf');
if ~isfolder(vox_results_dir)
    mkdir(vox_results_dir)
end
pethdr = spm_vol(vox_petfile); pethdr = pethdr(1);
vox_img_names = cell(length(vox_par_est_names), 1);

for i = 1:length(vox_par_est_names)
    if ~(isnan(model_options.filter_size) || model_options.filter_size == 0)
        pethdr.fname   = [vox_results_dir filesep subject '_' lower(char(vox_par_est_names(i))) '_' model_options.VOXmodel '_from_smoothed.nii'];
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
    temp_param_img(mask_idx) = vox_par_est(:,i);
    
    param(i,1)=mean(temp_param_img(cortexMaskFull==1));
    param(i,2)=mean(temp_param_img(cortexMaskFull==2));
    spm_write_vol(pethdr, temp_param_img);
    if contains(vox_par_est_names(i),'K1')
        vol=imgaussfilt3(temp_param_img,1);
        figure;
        tiledlayout(1,3,'TileSpacing','compact','Padding','tight')
        nexttile
        m=rot90(squeeze(max(cortexMaskFull,[],2)));
        I=rot90(squeeze(max(vol,[],2)));
        Ishow = mat2gray(I, [0 500]); 
        % Choose colors for labels 1 and 2
        cmap = [1 0 0;   % label 1 = red
                0 1 0];  % label 2 = green
        J = labeloverlay(Ishow, m, 'Colormap', cmap, 'Transparency', 0.45);
        imshow(J);
        title('MiP')
        nexttile
        [idxX, idxY, idxZ] = ind2sub(size(cortexMaskFull), find(cortexMaskFull));
        com_voxel = [mean(idxX), mean(idxY), mean(idxZ)];
        m=rot90(squeeze(cortexMaskFull(:,round(com_voxel(2)),:)));
        I=rot90(squeeze(vol(:,round(com_voxel(2)),:)));
        Ishow = mat2gray(I, [0 500]); 
        % Choose colors for labels 1 and 2
        cmap = [1 0 0;   % label 1 = red
                0 1 0];  % label 2 = green
        J = labeloverlay(Ishow, m, 'Colormap', cmap, 'Transparency', 0.45);
        imshow(J);
        title('Coronal')
        nexttile
        [idxX, idxY, idxZ] = ind2sub(size(cortexMaskFull), find(cortexMaskFull));
        com_voxel = [mean(idxX), mean(idxY), mean(idxZ)];
        m=rot90(squeeze(cortexMaskFull(:,:,round(com_voxel(3)))));
        I=rot90(squeeze(vol(:,:,round(com_voxel(3)))));
        Ishow = mat2gray(I, [0 500]); 
        % Choose colors for labels 1 and 2
        cmap = [1 0 0;   % label 1 = red
                0 1 0];  % label 2 = green
        J = labeloverlay(Ishow, m, 'Colormap', cmap, 'Transparency', 0.45);
        imshow(J);
        title('Axial')
        outdir = fullfile(subject_dir,'results','results_no_mc');
        saveas(gcf, fullfile(outdir, [subject '_kidney_VOIs.png']));
    end
end


%% ==== Save results ====
outdir = fullfile(subject_dir,'results','results_no_mc','voi-modelling','h2o');
if ~exist(outdir, 'dir'); mkdir(outdir); end

saveas(gcf, fullfile(outdir, [subject '_kidney_ROIfit.png']));
outdir = fullfile(subject_dir,'results','results_no_mc');
% You can rename columns later; for speed / robustness we keep generic names
Th = array2table([resultsH param(:,1)']);
Th.Properties.VariableNames={'K1 (voi)', 'k2 (voi)', 'Va (voi)', 'R2 (voi)', 'delay (voi)',...
    [vox_par_est_names{1} ' (vox)'],[vox_par_est_names{2} ' (vox)'],[vox_par_est_names{3} ' (vox)'],...
    [vox_par_est_names{4} ' (vox)'],[vox_par_est_names{5} ' (vox)'],[vox_par_est_names{6} ' (vox)']};
Tv = array2table([resultsV param(:,2)']);
Tv.Properties.VariableNames={'K1 (voi)', 'k2 (voi)', 'Va (voi)', 'R2 (voi)', 'delay (voi)',...
    [vox_par_est_names{1} ' (vox)'],[vox_par_est_names{2} ' (vox)'],[vox_par_est_names{3} ' (vox)'],...
    [vox_par_est_names{4} ' (vox)'],[vox_par_est_names{5} ' (vox)'],[vox_par_est_names{6} ' (vox)']};

writetable(Tv, fullfile(outdir,[subject '_kidney_results.xlsx']), "Sheet","Left");
writetable(Th, fullfile(outdir,[subject '_kidney_results.xlsx']), "Sheet","Right");

end
function M = refineAortaMask(M, petSum, se)
% Refine aorta mask using high-uptake threshold inside provided mask

if ~any(M(:))
    return;
end

M = imdilate(M, se);

vals = petSum(M);
th   = prctile(vals, 75);

M = (petSum .* M) > th;

M = imclose(M, se);
M = imfill(M, 'holes');
M = imerode(M, se);
M = logical(M);

end

function [x1,y1,z1, x2,y2,z2] = findBoundingBox(img)

img = greatestConnectedComponentOnly(img);
stats = regionprops(img, 'BoundingBox');
bb = stats(1).BoundingBox;   % [x y z dx dy dz] in voxel indices

x1 = round(bb(1));
y1 = round(bb(2));
z1 = round(bb(3));
x2 = x1 + round(bb(4)) - 1;
y2 = y1 + round(bb(5)) - 1;
z2 = z1 + round(bb(6)) - 1;

% Expand bounding box by 3 voxels each side (you can replace with physical margin)
margin = 3;
x1 = max(1, x1 - margin);
y1 = max(1, y1 - margin);
z1 = max(1, z1 - margin);

x2 = min(size(img,2), x2 + margin);
y2 = min(size(img,1), y2 + margin);
z2 = min(size(img,3), z2 + margin);

end

function largestComponent = greatestConnectedComponentOnly(img)

binaryVolume = img > 0;
cc = bwconncomp(binaryVolume, 6);

if cc.NumObjects == 0
    largestComponent = false(size(binaryVolume));
    return;
end

componentSizes = cellfun(@numel, cc.PixelIdxList);
[~, largestIdx] = max(componentSizes);

largestComponent = false(size(binaryVolume));
largestComponent(cc.PixelIdxList{largestIdx}) = true;

end

function study_specs=get_studyspecs(pet_orig_file)

json_file = strrep(pet_orig_file,'.nii','.json');

if exist(json_file,'file')
        json_str = fileread(json_file);
        % Convert JSON string to MATLAB variables
        json_metadata = jsondecode(json_str);
        frames = [];
        if isfield(json_metadata,'FrameTimesStart')
            frames(:,1) = json_metadata.FrameTimesStart - json_metadata.FrameTimesStart(1);         
        else 
             frames(:,1) = cumsum(json_metadata.FrameDuration) - json_metadata.FrameDuration;             
        end
        frames(:,2) = frames(:,1) + json_metadata.FrameDuration;
else
    error('json-file not found. Possible error during dicom to nifti conversion.')
end

study_specs = struct('dose', json_metadata.InjectedRadioactivity/1e6, 'scanner', json_metadata.ManufacturersModelName, ...
    'tracer', json_metadata.Radiopharmaceutical, 'weight', 0, 'frames', frames, 'glucose', 0, 'hct', NaN);

%assume all data is decay-corrected
study_specs.dc = 1;
end


function save_pmod_bld(TAC, filename)
    % TAC: N×2 matrix [time_in_seconds, activity_kBq_per_mL]
    % filename: e.g. 'mycurve.bld'
    
    % Extract columns
    time = TAC(:,1);
    activity = TAC(:,2);
    
    % Open file
    fid = fopen(filename, 'w');
    if fid < 0
        error('Cannot open file %s for writing.', filename);
    end
    
    % Write header (PMOD CRV / BLD format)
    fprintf(fid, 'time[seconds]\tblood[kBq/cc]\n');
    
    % Write each row
    for i = 1:length(time)
        fprintf(fid, '%.6f\t%.6f\n', time(i), activity(i));
    end
    
    fclose(fid);
    fprintf('Saved PMOD-compatible BLD file: %s\n', filename);
end
