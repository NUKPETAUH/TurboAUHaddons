function turbo_kidney_clean(folder, subject)
% TURBO_KIDNEY
%   Perform kidney ROI extraction and kinetic modeling for integration
%   in Turku PET Centre's TURBO pipeline using [15O]H2O PET with 
%   abdominal IDIFs.
%
%   turbo_kidney(folder, subject)
%
%   This function performs the following:
%       • Checks for .nii.gz files, unzips them only if necessary,
%         and cleans up temporary .nii files after processing.
%       • Loads PET volumes (raw + smoothed) and CT-based segmentations.
%       • Automatically extracts kidney bounding boxes and registers
%         kidney masks to PET space (affine imregtform).
%       • Generates renal cortical masks using morphology + PET uptake.
%       • Extracts TACs from each kidney cortex.
%       • Constructs 3 abdominal aorta IDIFs from multiple axial ranges.
%       • Reads an external cardiac IDIF and combines it with the
%         abdominal IDIFs.
%       • Runs ROI-level kinetic modeling using the standard
%         one-tissue water model (fitdelay + fit_h2o).
%       • Produces diagnostic plots of TAC fits.
%       • Writes Excel tables with K1, k2, Va, delay, and fit quality.
%       • Saves a PNG summary figure and displays a 3D aorta overlay.
%
%   INPUTS:
%       folder  : base directory of subject data (string)
%       subject : subject ID or name (string)
%
%   OUTPUTS:
%       Saves figures and Excel files to:
%           <folder>/nifti/<subject>/results/results_no_mc/vox-modelling/imgh2obf/
%
%   DEPENDENCIES:
%       fitdelay.m
%       fit_h2o.m
%       Image Processing Toolbox
%       Medical Imaging Toolbox (for medicalVolume)
%
%   NOTES:
%       • PET voxel spacing is detected automatically.
%       • All temporary .nii files extracted from .nii.gz are deleted
%         at the end to keep the folder clean.
%       • The function does not return data to the workspace; all
%         outputs are written to disk.
%
%   Author: Lars Tolbod
%   Date:   Nov 2025
% -------------------------------------------------------------------------
%% ==== Settings ====
showDebug = false;   % set true to show per-kidney volshow

%% ==== Close all figures ====
close all force;

%% === Unzip NIfTI files (only if needed) ===
disp('Checking NIfTI files...')

baseNifti = fullfile(folder, 'nifti', subject);
gzFiles = dir(fullfile(baseNifti, '**', '*.nii.gz'));
createdNiiFiles = {};   % track unzipped files created during this run

for k = 1:numel(gzFiles)
    gzPath = fullfile(gzFiles(k).folder, gzFiles(k).name);
    niiPath = gzPath(1:end-3);    % remove .gz → .nii

    needUnzip = false;

    if ~exist(niiPath, 'file')
        % No .nii → need to unzip
        needUnzip = true;
    else
        % Check timestamps: unzip only if .gz is newer than .nii
        gzInfo  = dir(gzPath);
        niiInfo = dir(niiPath);
        if gzInfo.datenum > niiInfo.datenum
            needUnzip = true;
        end
    end

    if needUnzip
        fprintf('Unzipping %s...\n', gzFiles(k).name);
        gunzip(gzPath);
        createdNiiFiles{end+1} = niiPath; %#ok<AGROW>
    end
end


%% ==== Load data ====
disp('Loading data...');
segV       = medicalVolume(fullfile(baseNifti, 'anat', [subject '_ct_segmentations_coregistered.nii']));
petV       = medicalVolume(fullfile(baseNifti, 'pet',  [subject '_pet.nii']));

if isfile(fullfile(baseNifti, 'pet',  [subject '_pet_smoothed.nii']))

    petVsmooth = medicalVolume(fullfile(baseNifti, 'pet',  [subject '_pet_smoothed.nii']));
else
    disp('Smoothed 4D not found - applying gaussian filter')
    petVsmooth=petV;
    V = petVsmooth.Voxels;                  % voxel data
    sp = petVsmooth.VoxelSpacing;           % [sx sy sz] in mm
    target_mm = 3;                    % desired Gaussian sigma in mm
    
    % sigma per spatial dimension in voxels
    sigma_vox = target_mm ./ sp;      % 1x3 vector
    
    % allocate filtered volume with same class
    Vf = zeros(size(V), 'like', V);
    
    num_frames = size(V,4);
    for frame = 1:num_frames
        Vf(:,:,:,frame) = imgaussfilt3(V(:,:,:,frame), sigma_vox);
    end
    petVsmooth.Voxels=Vf;
    clear V sp target_mm num_frames frame
end

turbo_options = load(fullfile(baseNifti, 'options', [subject '_turbo_options.mat']));
model_options = load(fullfile(baseNifti, 'options', [subject '_modelling_options.mat']));
study_options = load(fullfile(baseNifti, 'options', [subject '_study_specs.mat']));

roiNames = turbo_options.roinames;

% ---- Automatic voxel spacing detection ----
% medicalVolume.VoxelSpacing is [dy dx dz] or similar; we only need
% the physical size in each dimension and sometimes the mean.
voxelSpacing = petV.VoxelSpacing(:)';  % [sY sX sZ] (mm)
voxelSizeMean = mean(voxelSpacing);    % mm, used for isotropic morphology
dz = voxelSpacing(3);                  % slice thickness in mm

%% ==== Identify ROIs ====
containsKidney = contains(roiNames, "kidney", 'IgnoreCase', true);
containsCyst   = contains(roiNames, "cyst",   'IgnoreCase', true);
kidney_label   = find(containsKidney & ~containsCyst);
aorta_label    = find(contains(roiNames, "aorta", 'IgnoreCase', true));

numKidneys = numel(kidney_label);
if numKidneys == 0
    error('No kidney labels found in turbo_options.roinames.');
end

%% ==== Registration optimizer ====
[optimizer, metric] = imregconfig('multimodal');
optimizer.InitialRadius       = 0.001;
optimizer.GrowthFactor        = 1.005;
optimizer.MaximumIterations   = 400;

%% ==== Pre-extract PET data once (speed) ====
petData      = petV.Voxels;        % 4D, large but already in memory anyway
petDataSmooth = petVsmooth.Voxels; % 4D smoothed

kidneyPET    = cell(1, numKidneys);
cortexMasks  = cell(1, numKidneys);
zBounds      = zeros(2, numKidneys);

%% ==== Process each kidney ====
disp('Extracting cortical masks and kidney TACs...');
for ki = 1:numKidneys

    roiName = roiNames{kidney_label(ki)};
    isRightKidney = contains(roiName, 'right', 'IgnoreCase', true);

    % Kidney binary mask from segmentation
    voxelMask = (segV.Voxels == kidney_label(ki));

    % Bounding box in index space
    [x1,y1,z1, x2,y2,z2] = findBoundingBox(voxelMask);
    zBounds(:,ki) = [z1; z2];

    % Crop smoothed PET & kidney mask
    petCrop     = petDataSmooth(y1:y2, x1:x2, z1:z2, :);
    petMeanCrop = mean(petCrop, 4);
    maskCrop    = voxelMask(y1:y2, x1:x2, z1:z2);

    % Registration: kidney mask -> PET
    maskSmooth = imgaussfilt3(single(maskCrop), 2.0);
    tform = imregtform(maskSmooth, petMeanCrop, 'affine', optimizer, metric);
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
    kidneyPET{ki}   = petData(y1:y2, x1:x2, z1:z2, :);
    cortexMasks{ki} = finalMask;

    % Optional visualization per kidney
    if showDebug
        viewer2 = viewer3d;
        volshow(petMeanCrop, 'OverlayData', finalMask, ...
            'OverlayRenderingStyle','LabelOverlay', ...
            'RenderingStyle','SlicePlanes', ...
            'OverlayAlpha',0.7, 'Parent',viewer2);
        title(viewer2, roiName);
    end

    % You currently compute hilum stuff elsewhere; spacing is ready:
    %   voxelSpacing, isRightKidney are available here if needed.
    %#ok<NASGU> % isRightKidney not used yet

end

%% ==== Aorta mask generation (using physical dz) ====
disp('Building aorta masks...');
petSum = sum(petData, 4);     % 3D sum over time

kidneyZmin = min(zBounds(1,:));
kidneyZmax = max(zBounds(2,:));
kidneySpanSlices = kidneyZmax - kidneyZmin + 1;
kidneyDepth_mm   = kidneySpanSlices * dz;

targetSAFOV_mm = 250; % 25 cm
padSlices = max(0, round((targetSAFOV_mm - kidneyDepth_mm) / (2*dz)));
SAFOV = [max(1, kidneyZmin - padSlices), ...
         min(size(segV.Voxels,3), kidneyZmax + padSlices)];

aortaMask = cell(1,4);

baseMask = false(size(segV.Voxels));
baseMask(segV.Voxels == aorta_label) = true;
seAorta = strel('cube',3);

% Mask 1 = original
aortaMask{1} = baseMask;

% Mask 2: full SAFOV slice range
aortaMask{2} = false(size(baseMask));
aortaMask{2}(:,:,SAFOV(1):SAFOV(2)) = baseMask(:,:,SAFOV(1):SAFOV(2));
aortaMask{2} = refineAortaMask(aortaMask{2}, petSum, seAorta);

% Mask 3: upper half SAFOV
midZ = round(mean(SAFOV));
aortaMask{3} = false(size(baseMask));
aortaMask{3}(:,:,midZ:SAFOV(2)) = baseMask(:,:,midZ:SAFOV(2));
aortaMask{3} = refineAortaMask(aortaMask{3}, petSum, seAorta);

% Mask 4: above kidneys
aortaMask{4} = false(size(baseMask));
aortaMask{4}(:,:,kidneyZmax:SAFOV(2)) = baseMask(:,:,kidneyZmax:SAFOV(2));
aortaMask{4} = refineAortaMask(aortaMask{4}, petSum, seAorta);

%% ==== TAC extraction (kidneys) ====
disp('Extracting kidney TACs...');
numFrames = size(petData,4);
tac = zeros(numFrames, numKidneys);

for ki = 1:numKidneys
    mask = cortexMasks{ki};
    maskIdx = mask(:);   % speed: linear indexing

    % reshape kidneyPET to 2D [nVox x nFrames] for vectorized mean
    vol4D = kidneyPET{ki};
    vol2D = reshape(vol4D, [], size(vol4D,4));   % [nVox x nFrames]
    tac(:,ki) = mean(vol2D(maskIdx,:), 1).';     % mean over voxels
end

%% ==== IDIF extraction (aorta) ====
disp('Extracting IDIFs...');
idif = zeros(numFrames, 3);
for f = 1:numFrames
    frame = petData(:,:,:,f);
    for k = 1:3
        idif(f,k) = mean(frame(aortaMask{k+1}));
    end
end
inp_file=fullfile(baseNifti, 'input', [subject '_IDIFaorta_no_mc.txt']);
if ~isfile(inp_file)
    inp_file=fullfile(baseNifti, 'input', [subject '_IDIFaorta_mc_rigid.txt']);
    if ~isfile(inp_file)
        disp('Input file not found')
        
    end
end
T = readtable(inp_file);
idif = [T.Var2, idif];  % first column TURBO cardiac IDIF
idif_label = {'TURBO Cardiac IDIF','IDIF1','IDIF2','IDIF3'};

%% ==== ROI-level modeling ====
disp('Performing ROI-level modeling...');

idif = idif ./ 1000;
tac  = tac  ./ 1000;

figure;
tiledlayout(3,4);

% IDIF overview
ax1 = nexttile([1 2]);
plot(ax1, T.Var1, idif);
legend(idif_label, 'Location','best');
title(ax1, 'IDIFs');
xlabel(ax1,'Time (min)');
ylabel(ax1,'Activity (kBq/mL)');

ax2 = nexttile(3,[1 2]);
plot(ax2, T.Var1, idif);
legend(ax2, idif_label, 'Location','best');
xlim(ax2, [0 1.5]);
title(ax2, 'Early IDIFs');
xlabel(ax2,'Time (min)');
ylabel(ax2,'Activity (kBq/mL)');

frametimes = study_options.frames * 60;
inputTimes = mean(study_options.frames,2) * 60;

nIF    = size(idif,2);
nParam = 5;      % [K1 k2 Va R2 delay]
resultsH = zeros(1, nIF*nParam);
resultsV = zeros(1, nIF*nParam);

for id = 1:nIF
    tacIF = idif(:,id);
    inputdata = double([inputTimes tacIF]);

    for ki = 1:numKidneys

        [~, corrected, fitdelay_result, ~, ~, ~] = fitdelay( ...
            frametimes, tac(:,ki)', inputdata, 1, 'verbose', 5);

        [modelfit, x_opt, ~, ~, ~] = fit_h2o( ...
            frametimes, tac(:,ki), ...
            corrected(:,1), corrected(:,2), corrected(:,2), ...
            0, 1, ...
            model_options.n_iter, ...
            model_options.lb(1,:), model_options.ub(1,:));

        % Extract parameters
        K1   = x_opt(1) * 100 * 60;
        k2   = (x_opt(1)/x_opt(2)) * 60;
        Va   = x_opt(3) * 100;
        dly  = fitdelay_result + x_opt(4);
        flow = (1 - x_opt(3)) * x_opt(1) * 100 * 60; %#ok<NASGU>

        r  = corrcoef(tac(:,ki), modelfit);
        R2 = r(1,2)^2;

        % Plot model fit
        nexttile;
        plot(inputdata(:,1), modelfit, 'b', ...
             inputdata(:,1), tac(:,ki), 'x', ...
             inputdata(:,1), inputdata(:,2)*Va/100, 'r');
        title({[strrep(roiNames{kidney_label(ki)},'_',' ') '  ' idif_label{id}], ...
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

%% ==== Save results ====
outdir = fullfile(baseNifti,'results','results_no_mc','roi-modelling','fit_h2o','kidney');
if ~exist(outdir, 'dir'); mkdir(outdir); end

saveas(gcf, fullfile(outdir, [subject '_kidney_ROIfit.png']));

% You can rename columns later; for speed / robustness we keep generic names
Th = array2table(resultsH);
Tv = array2table(resultsV);

writetable(Tv, fullfile(outdir,[subject '_kidney_ROIfit.xlsx']), "Sheet","Left");
writetable(Th, fullfile(outdir,[subject '_kidney_ROIfit.xlsx']), "Sheet","Right");

%% === Clean up created NIfTI files ===
if exist('createdNiiFiles', 'var') && ~isempty(createdNiiFiles)
    disp('Cleaning up temporary .nii files...');
    for i = 1:numel(createdNiiFiles)
        if exist(createdNiiFiles{i}, 'file')
            delete(createdNiiFiles{i});
        end
    end
end

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
