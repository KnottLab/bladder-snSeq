clear all
close all force
addpath('~/Documents/MATLAB/','~/Documents/MATLAB/functions/cbrewer/')

samples = {'B0','B1','B2','B3','B4','B5A','B5B','B6A','B6B','B7A','B7B','B8A','B8B'};

for i = 1:size(samples,2)

    metadata = readtable(['~/Documents/Bladder/seuratOutput/',samples{i},'_preDemux_metadata.csv'],'ReadRowNames',true,'ReadVariableNames',true);
    dataMat = readtable(['~/Documents/Bladder/seuratOutput/',samples{i},'_preDemux_HASH_data.csv'],'ReadRowNames',true,'ReadVariableNames',true);

    YourArray = table2array(dataMat);
    YourNewTable = array2table(YourArray.');
    YourNewTable.Properties.RowNames = dataMat.Properties.VariableNames;
    YourNewTable.Properties.VariableNames = matlab.lang.makeValidName(dataMat.Properties.RowNames);
    combinedDataMat = join(metadata,YourNewTable,'Keys','RowNames');

    %recovering negative cells
    
    neg = combinedDataMat(strcmp(combinedDataMat.HASH_classification_global,'Negative'),:);
    sing = combinedDataMat(strcmp(combinedDataMat.HASH_classification_global,'Singlet'),:);

    if strcmp(samples{i},'B0')
        HASH1cutoff = min(sing.Hash1_TTCCTGCCATTACTA(strcmp(sing.HASH_maxID,'Hash1-TTCCTGCCATTACTA')));
        HASH2cutoff = min(sing.Hash2_CCGTACCTCATTGTT(strcmp(sing.HASH_maxID,'Hash2-CCGTACCTCATTGTT')));
        HASH3cutoff = min(sing.Hash3_GGTAGATGTCCTCAG(strcmp(sing.HASH_maxID,'Hash3-GGTAGATGTCCTCAG')));
    
        neg.HASH1pos = zeros(size(neg,1),1);
        neg.HASH1pos(neg.Hash1_TTCCTGCCATTACTA>HASH1cutoff) = 1;
        neg.HASH2pos = zeros(size(neg,1),1);
        neg.HASH2pos(neg.Hash2_CCGTACCTCATTGTT>HASH2cutoff) = 1;
        neg.HASH3pos = zeros(size(neg,1),1);
        neg.HASH3pos(neg.Hash3_GGTAGATGTCCTCAG>HASH3cutoff) = 1;
        
        neg.totalPos = sum(table2array(neg(:,17:end)),2);
        recoveredNegs = neg(neg.totalPos==1,1:end-7);
    else
        HASH1cutoff = min(sing.Hash1_TTCCTGCCATTACTA(strcmp(sing.HASH_maxID,'Hash1-TTCCTGCCATTACTA')));
        HASH2cutoff = min(sing.Hash2_CCGTACCTCATTGTT(strcmp(sing.HASH_maxID,'Hash2-CCGTACCTCATTGTT')));
        HASH3cutoff = min(sing.Hash3_GGTAGATGTCCTCAG(strcmp(sing.HASH_maxID,'Hash3-GGTAGATGTCCTCAG')));
        HASH4cutoff = min(sing.Hash4_TGGTGTCATTCTTGA(strcmp(sing.HASH_maxID,'Hash4-TGGTGTCATTCTTGA')));
        HASH5cutoff = min(sing.Hash5_ATGATGAACAGCCAG(strcmp(sing.HASH_maxID,'Hash5-ATGATGAACAGCCAG')));
        HASH6cutoff = min(sing.Hash6_CTCGAACGCTTATCG(strcmp(sing.HASH_maxID,'Hash6-CTCGAACGCTTATCG')));

        neg.HASH1pos = zeros(size(neg,1),1);
        neg.HASH1pos(neg.Hash1_TTCCTGCCATTACTA>HASH1cutoff) = 1;
        neg.HASH2pos = zeros(size(neg,1),1);
        neg.HASH2pos(neg.Hash2_CCGTACCTCATTGTT>HASH2cutoff) = 1;
        neg.HASH3pos = zeros(size(neg,1),1);
        neg.HASH3pos(neg.Hash3_GGTAGATGTCCTCAG>HASH3cutoff) = 1;
        neg.HASH4pos = zeros(size(neg,1),1);
        neg.HASH4pos(neg.Hash4_TGGTGTCATTCTTGA>HASH4cutoff) = 1;
        neg.HASH5pos = zeros(size(neg,1),1);
        neg.HASH5pos(neg.Hash5_ATGATGAACAGCCAG>HASH5cutoff) = 1;
        neg.HASH6pos = zeros(size(neg,1),1);
        neg.HASH6pos(neg.Hash6_CTCGAACGCTTATCG>HASH6cutoff) = 1;
        
        neg.totalPos = sum(table2array(neg(:,20:end)),2);
        recoveredNegs = neg(neg.totalPos==1,1:end-13);
    end

    sum(neg.totalPos==1)
    
    %fix metadata for recovered cells
    
    if i==1
        filteredmetadata = readtable('~/Documents/Bladder/seuratOutput/combinedSobj_filtered_metadata.csv','ReadRowNames',true,'ReadVariableNames',true);
    else
        filteredmetadata = readtable('~/Documents/Bladder/seuratOutput/combinedSobj_filtered_metadata_plusRecovered.csv','ReadRowNames',true,'ReadVariableNames',true);
    end
    
    recoveredNegs.Patient_Rep = repmat({''},size(recoveredNegs,1),1);
    recoveredNegs.Patient = repmat({''},size(recoveredNegs,1),1);
    recoveredNegs.cohort = repmat({''},size(recoveredNegs,1),1);

    for j = 1:size(recoveredNegs,1)
        match = filteredmetadata(strcmp(filteredmetadata.sampleID,recoveredNegs.sampleID{i})&strcmp(filteredmetadata.HASH_maxID,recoveredNegs.HASH_maxID{i}),:);
        if size(match,1) == 0
            recoveredNegs.Patient_Rep(j) = {'NA'};
            test = strsplit(recoveredNegs.HASH_maxID{j},'-');
            recoveredNegs.Patient(j) = test(1);
            recoveredNegs.cohort(j) = {'NEW'};
        else
            recoveredNegs.Patient_Rep(j) = match.Patient_Rep(1);
            recoveredNegs.Patient(j) = match.Patient(1);
            recoveredNegs.cohort(j) = match.cohort(1);
        end
    end

    recoveredNegs.Properties.RowNames = strcat(recoveredNegs.sampleID,'__',recoveredNegs.Properties.RowNames);

    filteredmetadata = [filteredmetadata;recoveredNegs];
    writetable(filteredmetadata,'~/Documents/Bladder/seuratOutput/combinedSobj_filtered_metadata_plusRecovered.csv','WriteRowNames',true,'WriteVariableNames',true)
    
end
