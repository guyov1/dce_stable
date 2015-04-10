function [Mask,data4D_filtered,Filter,Filter_total] = data4D_mask_filtering(voxels_data_4D,Mask,filter_is_on)

% MASK1 - thresholding to delete blank voxels:
% We separate voxels to "data" (over th) and "noise" (below th)
% Based on the first time sample.

[N_rows,N_cols,N_slices,N_time_points]=size(voxels_data_4D);

%% loop on the filters:
filter_names={'LowAtZero';'LowAverageBaseline';'WrongBolusPeak';'PeakSaturation';'BigFluctuations';'ZeroValues';'LowSteadyState'};
% filter_is_on=[1 1 0 1 1];

N_filters=length(filter_names);

filter_thresholds=[500 500 0 0.5 1.1 1 0];

filter_show_stats=[0 0 0 0 0 1 0];
filter_show_mask=[0 0 0 0 0 0 0];

Filter=struct;
Filter(1).mean_DSC_curve=mean(voxels_data_4D,4);
Filter(1).first_time_sample_to_use_next=5;
Filter(2).first_baseline_point_to_include_total=10;
Filter(2).first_baseline_point_to_include=Filter(2).first_baseline_point_to_include_total-Filter(1).first_time_sample_to_use_next+1;
Filter(2).baseline_length_for_mask=10;
Filter(3).first_baseline_point_to_include_total=10;
Filter(3).first_baseline_point_to_include=Filter(3).first_baseline_point_to_include_total-Filter(1).first_time_sample_to_use_next+1;
Filter(3).baseline_length_for_mask=10;
Filter(3).bolus_range_from_mean=6;
Filter(3).bolus_range_length=40;
Filter(4).bolus_range_from_mean=8;
Filter(7).first_baseline_point_to_include_total=10;
Filter(7).first_baseline_point_to_include=Filter(7).first_baseline_point_to_include_total-Filter(1).first_time_sample_to_use_next+1;
Filter(7).baseline_length_for_mask=10;
Filter(7).last_point_to_include=size(voxels_data_4D,4)-5;
Filter(7).steady_state_part_length=25;



N_time_points_init=N_time_points;
voxels_data_4D_curr=voxels_data_4D;
% initializing filter data, in case no filter is used at all
data_for_filter=ones(N_rows,N_cols,N_slices);
th=0;
Mask_prev=Mask;
[data_voxels_row_col_slice,data_voxels_indices,noise_voxels_indices,Mask] = filter_thresholding(data_for_filter,th,Mask_prev,[]);
Filter_total.data_voxels_row_col_slice=data_voxels_row_col_slice;
Filter_total.data_voxels_indices=data_voxels_indices;

Mask_prev=Mask;
for filt_ind=1:N_filters
    if filter_is_on(filt_ind)
        Filter(filt_ind).name=filter_names{filt_ind};
        display(['     ',Filter(filt_ind).name]);
        Filter(filt_ind).th=filter_thresholds(filt_ind);
        
        %%%%%% Preprocess to create a specific data needed for the filtering:
        filter_str=['PreProcess_',Filter(filt_ind).name];
        filter_preprocess_func=str2func(filter_str);
        [data_for_filter ref] = filter_preprocess_func(voxels_data_4D_curr,Filter(filt_ind),Mask_prev);
        
        
        %%%%%% Thresholding to create new mask: (we use the data after pre-process)
        th=filter_thresholds(filt_ind);
        [data_voxels_row_col_slice,data_voxels_indices,noise_voxels_indices,Mask] = filter_thresholding(data_for_filter,th,Mask_prev,ref);
        Filter(filt_ind).data_voxels_row_col_slice=data_voxels_row_col_slice;
        Filter(filt_ind).data_voxels_indices=data_voxels_indices;
        N_data_voxels=length(data_voxels_indices);
        Filter_total.data_voxels_row_col_slice=data_voxels_row_col_slice;
        Filter_total.data_voxels_indices=data_voxels_indices;
        %%%%%% Applying the new mask on the 4D data:
        voxels_data_4D_next=zeros(size(voxels_data_4D_curr));
        for t=1:N_time_points
            voxels_data_4D_next(:,:,:,t)=voxels_data_4D_curr(:,:,:,t).*Mask;
        end
        
        
        %%%%%% Statistics for the data voxels after applying the new Mask:
        mean_DSC_curve=zeros(1,N_time_points_init);
        std_DSC_curve=zeros(1,N_time_points_init);
        for t=1:N_time_points_init
            voxels_data_3D_mask_t=voxels_data_4D(:,:,:,t);
            voxels_data_3D_mask_t_vec=voxels_data_3D_mask_t(Mask>0);
            mean_DSC_curve(t)=mean(voxels_data_3D_mask_t_vec);
            std_DSC_curve(t)=std(voxels_data_3D_mask_t_vec);
        end
        
        Filter(filt_ind+1).mean_DSC_curve=mean_DSC_curve;
        %plot some random DSC curves and the mean of ALL DSC curves
        N_DSC_curves_to_show=20;
        check_vec=sort(random('unid',N_data_voxels,1,N_DSC_curves_to_show));
        data_voxels_rows=data_voxels_row_col_slice(:,1);
        data_voxels_cols=data_voxels_row_col_slice(:,2);
        data_voxels_slices=data_voxels_row_col_slice(:,3);
        
        if filter_show_stats(filt_ind)
            figure;
            for jj=1:length(check_vec)
                check_ind=check_vec(jj);
                plot(squeeze(voxels_data_4D(data_voxels_rows(check_ind),data_voxels_cols(check_ind),data_voxels_slices(check_ind),:)));
                hold all;
            end
            std_part_to_show=1;
            title([num2str(N_DSC_curves_to_show),' random DSC curves out of ',num2str(N_data_voxels),' voxels passed "',filter_names{filt_ind},'" filter. blue bold = mean of passed voxels + ',num2str(std_part_to_show),'std. Average std (over time) =',num2str(mean(std_DSC_curve)),'. TH=',num2str(th)]);
            %             errorbar(mean_DSC_curve,std_part_to_show*std_DSC_curve,'--b','LineWidth',2);
            plot(mean_DSC_curve,'--b','LineWidth',3);
            %     plot(mean_DSC_curve_total,'-.g','LineWidth',3);
            hold off;
        end
        
        % Show the mask on the brain images 
        if filter_show_mask(filt_ind)
%             voxels_data_to_show=voxels_data_4D_next(:,:,:,);
%             voxels_data_to_show=repmat(voxels_data_4D_next(:,:,:,floor(N_time_points/2)),[1,1,1,3]);
%             voxels_data_to_show(noise_voxels_indices)=[1 0 0]; %red color
%             voxels_data_to_show(data_voxels_indices)=[1 1 1]*voxels_data_to_show(data_voxels_indices);
            figure(500);
            title(['Mask after the filter ',Filter(filt_ind).name]);
            for slice=1:N_slices
               subplot(floor(N_slices/4),5,slice);
               imshow(mat2gray(voxels_data_4D_next(:,:,slice,floor(N_time_points/2))));
            end
        end
    else  % filter is off
        Mask=Mask_prev;
        voxels_data_4D_next=voxels_data_4D_curr;
        Filter(filt_ind+1).mean_DSC_curve=Filter(filt_ind).mean_DSC_curve;
    end
    
    % truncate the first time samples, after applying the first filter:
    if filt_ind==1
        voxels_data_4D_next=voxels_data_4D_curr(:,:,:,Filter(filt_ind).first_time_sample_to_use_next:end);
        N_time_points=N_time_points-Filter(filt_ind).first_time_sample_to_use_next+1;
    end
    
    %%%%%% Substitution before next iteration of the loop:
    Mask_prev=Mask;
    voxels_data_4D_curr=voxels_data_4D_next;
    
end

display('Done calculating filters.');

% arrange the output results:
% data4D_filtered=voxels_data_4D_curr;
data4D_filtered=voxels_data_4D.*repmat(Mask,[1 1 1 size(voxels_data_4D,4)]);;
% Mask is already ready



end
