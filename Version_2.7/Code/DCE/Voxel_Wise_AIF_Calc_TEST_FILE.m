
check=loadniidata('D:\users\guyn\DCE_OUT\ShTo_20081116\DCEMean.nii');


figure;
subplot(2,2,1);
imshow(mat2gray(DSC_Map_norm_final_rep_voxel(:,:,Slice_num_to_disp)));
subplot(2,2,3);
imshow(2*mat2gray(check(:,:,Slice_num_to_disp)));
subplot(2,2,2);
imshow(mat2gray(DSC_Map_norm_final_rep_voxel(:,:,Slice_num_to_disp+1)));
subplot(2,2,4);
imshow(2*mat2gray(check(:,:,Slice_num_to_disp+1)));


MTT_masked_2D = DSC_Map_norm_final_rep_voxel(DBrainMask);

regular_sigma_1 = OutAIFParam(3);
% Initiate AIF matrices
HAIF_local  = zeros(numel(MTT_masked_2D),size(HAIF,1),size(HAIF,2));
CHAIF_local = zeros(numel(MTT_masked_2D),size(HAIF,1),size(HAIF,2));
SAIF_local  = zeros([numel(MTT_masked_2D) nTDif numel(SampleTs)]);
CSAIF_local = zeros([numel(MTT_masked_2D) nTDif numel(SampleTs)]);

% Measure calculation time
tic;
for idx = 1:numel(MTT_masked_2D)
    
    % Measure time after 100 iterations and estimate finish time
    if (idx == 100)
        time_for_100_iter = toc;s
        display(sprintf('-I- Caclculating AIF for 100 voxels took %.2f seconds.',time_for_100_iter));
        display(sprintf('-I- Number of total voxels is %.2f.',numel(MTT_masked_2D)));
        time_for_total_min = ( numel(MTT_masked_2D) / 100 ) * time_for_100_iter / 60;
        display(sprintf('-I- Estimated time for voxel wise AIF calculation is %.2f minutes...',time_for_total));    
    end
    
    % Change sigma only if bigger than 0
    if ( MTT_masked_2D(idx) > 0 ) 
        modified_sigma_1 = MTT_masked_2D(idx) * regular_sigma_1;
        OutAIFParam(3)   = modified_sigma_1;
    else
        OutAIFParam(3)   = regular_sigma_1;
    end
    
    % Create Parker's AIF
    HAIF_local(idx,:,:) = AIF_Parker9t(OutAIFParam,HSampleTs);
    % Create Cummulative Parker's AIF
    CHAIF_local(idx,:,:) = cumtrapz(HSampleTs,HAIF_local(idx,:,:));
    
    
    for i=1:nTDif
        SAIF_local(idx,i,:)  = interp1(HSampleTs, squeeze(HAIF_local(idx,:,:)) ,SampleTs+TDif(i), [], 'extrap');
        CSAIF_local(idx,i,:) = interp1(HSampleTs, squeeze(CHAIF_local(idx,:,:)),SampleTs+TDif(i), [], 'extrap');
    end
    
end
