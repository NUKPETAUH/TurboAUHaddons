clear;clc
folder='/home/turbo/Documents/data/448/nifti/';
subjects=dir(folder);
pyexe='/home/turbo/turboSW/rhIDIF/rhIDIF_env/bin/python';
idifFolder='/home/turbo/turboSW/IDIF';


for i=3:length(subjects)
    disp(subjects(i).name)
    % Move original turbo brain folder
    brainFolder=[subjects(i).folder filesep subjects(i).name '/pet/pet_mc_rigid/brain'];
    if isfolder(brainFolder)
        disp('Moving original pet brain folder')
        movefile(brainFolder,[brainFolder '_org'])
    end
    magiaResultsFolder=[subjects(i).folder filesep subjects(i).name '/results/results_mc_rigid/magia_brain_results'];
    if isfolder(brainFolder)
        disp('Moving original results magia folder')
        movefile(bmagiaResultsFolder,[magiaResultsFolder '_org'])
    end
    inputFolder=[subjects(i).folder filesep subjects(i).name '/input'];
    if isfolder(inputfolder) && ~isfolder([inputfolder '/org'])
        disp('Moving original results magia folder')
        movefile(inputfolder,[inputfolder '/org'])
    end

    % Run RH IDIF python scipt
    petFile=[subjects(i).folder filesep subjects(i).name '/pet/' subjects(i).name '_pet.nii.gz'];
    segFile=[subjects(i).folder filesep subjects(i).name '/anat/' subjects(i).name '_ct_segmentations_coregistered.nii.gz'];
    outFolder=[subjects(i).folder filesep subjects(i).name '/input/rhIDIF'];
    if isfile(petFile) && isfile(segFile)
        command=[pyexe ' "' idifFolder filesep 'idif_rsl.py"' ' -i "' petFile '" -s "' segFile '" -o "' outFolder '"'];
        disp(command)
        [status,cmdout] = system(command,'-echo');
        if status>0
            break
        end
    else
        disp('Files for IDIF segmentation not found')
    end
  
    % Read RH IDIF 3 and 4
    t3=readtable([outFolder 'IDIF_totalsegmentator_segment-3.tac'],'FileType','text');
    t4=readtable([outFolder 'IDIF_totalsegmentator_segment-4.tac'],'FileType','text');      

    % average 3 and 4 and save as input file for turbo
    newT=array2table([(t3.start_seconds_+t3.end_Bq_cc_)/2/60 (t3.x__idif__+t4.x__idif__)/2]);
    outIDIFFile=[subjects(i).folder filesep subjects(i).name '/input/' subjects(i).name '_IDIFaorta_mc_rigid.txt'];
    writetable(newT,outIDIFFile,'WriteVariableNames',false,'Delimiter',',')
    disp(['New file written: ' outIDIFFile])

end