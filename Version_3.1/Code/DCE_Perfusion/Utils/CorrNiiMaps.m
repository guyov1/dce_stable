function [Rval,  Pval] = CorrNiiMaps( Path_Map1, Path_Map2, SliceNum, MinRange, MaxRange, ...
    Mask1, Mask2, Mask3, Mask4, Use_Mask, Mask_Names, Remove_Extreme_Vals, font_size, font_legend, map_name, OutPath, Close_CorrPlot_Fig, Remove_BV, Scatter_no_legend, Subj_Name)


NiiMap1   = loadniidata(Path_Map1);
NiiMap2   = loadniidata(Path_Map2);

NiiMap1   = NiiMap1(:, :, SliceNum);
NiiMap2   = NiiMap2(:, :, SliceNum);

Mask1 = uint8(Mask1~=0);
Mask2 = uint8(Mask2~=0);
Mask3 = uint8(Mask3~=0);
Mask4 = uint8(Mask4~=0);

Mask = Mask1 | Mask2 | Mask3 | Mask4;

NiiMap1_Masked = zeros(size(NiiMap1));
NiiMap2_Masked = zeros(size(NiiMap2));

NiiMap1_Masked1= zeros(size(NiiMap1));
NiiMap1_Masked2= zeros(size(NiiMap1));
NiiMap1_Masked3= zeros(size(NiiMap1));
NiiMap1_Masked4= zeros(size(NiiMap1));

NiiMap2_Masked1= zeros(size(NiiMap2));
NiiMap2_Masked2= zeros(size(NiiMap2));
NiiMap2_Masked3= zeros(size(NiiMap2));
NiiMap2_Masked4= zeros(size(NiiMap2));

if Use_Mask
    
    NiiMap1_Masked(Mask(:, :, SliceNum)==1)   = NiiMap1(Mask(:, :, SliceNum)==1);
    NiiMap2_Masked(Mask(:, :, SliceNum)==1)   = NiiMap2(Mask(:, :, SliceNum)==1);
    
    figure;
    subplot(2,2,1);
    imshow(NiiMap1);
    title(['DCE Map. Slice: ' num2str(SliceNum)]);
    subplot(2,2,2);
    imshow(NiiMap1_Masked,'Colormap',jet(255));
    title(['DCE Masked (All Segments) ']);
    subplot(2,2,3);
    imshow(NiiMap2);
    title(['DSC Map. Slice: ' num2str(SliceNum)]);
    subplot(2,2,4);
    imshow(NiiMap2_Masked,'Colormap',jet(255));
    title(['DSC Masked (All Segments) ']);
    
    % Per segment mask
    NiiMap1_Masked1(Mask1(:, :, SliceNum)==1)   = NiiMap1(Mask1(:, :, SliceNum)==1);
    NiiMap2_Masked1(Mask1(:, :, SliceNum)==1)   = NiiMap2(Mask1(:, :, SliceNum)==1);
    NiiMap1_Masked2(Mask2(:, :, SliceNum)==1)   = NiiMap1(Mask2(:, :, SliceNum)==1);
    NiiMap2_Masked2(Mask2(:, :, SliceNum)==1)   = NiiMap2(Mask2(:, :, SliceNum)==1);
    NiiMap1_Masked3(Mask3(:, :, SliceNum)==1)   = NiiMap1(Mask3(:, :, SliceNum)==1);
    NiiMap2_Masked3(Mask3(:, :, SliceNum)==1)   = NiiMap2(Mask3(:, :, SliceNum)==1);
    NiiMap1_Masked4(Mask4(:, :, SliceNum)==1)   = NiiMap1(Mask4(:, :, SliceNum)==1);
    NiiMap2_Masked4(Mask4(:, :, SliceNum)==1)   = NiiMap2(Mask4(:, :, SliceNum)==1);
    
    figure;
    subplot(2,4,1);
    imshow(NiiMap1_Masked1,'Colormap',jet(255));
    title(['DCE Masked (' char(Mask_Names{1}) ')']);
    subplot(2,4,2);
    imshow(NiiMap1_Masked2,'Colormap',jet(255));
    title(['DCE Masked (' char(Mask_Names{2}) ')']);
    subplot(2,4,3);
    imshow(NiiMap1_Masked3,'Colormap',jet(255));
    title(['DCE Masked (' char(Mask_Names{3}) ')']);
    subplot(2,4,4);
    imshow(NiiMap1_Masked4,'Colormap',jet(255));
    title(['DCE Masked (' char(Mask_Names{4}) ')']);
    subplot(2,4,5);
    imshow(NiiMap2_Masked1,'Colormap',jet(255));
    title(['DSC Masked (' char(Mask_Names{1}) ')']);
    subplot(2,4,6);
    imshow(NiiMap2_Masked2,'Colormap',jet(255));
    title(['DSC Masked (' char(Mask_Names{2}) ')']);
    subplot(2,4,7);
    imshow(NiiMap2_Masked3,'Colormap',jet(255));
    title(['DSC Masked (' char(Mask_Names{3}) ')']);
    subplot(2,4,8);
    imshow(NiiMap2_Masked4,'Colormap',jet(255));
    title(['DSC Masked (' char(Mask_Names{4}) ')']);
    
else
    NiiMap1_Masked = NiiMap1;
    NiiMap2_Masked = NiiMap2;
end

maxVal1   = max(max(max(NiiMap1_Masked)));
maxVal2   = max(max(max(NiiMap2_Masked)));
minVal1   = min(min(min(NiiMap1_Masked)));
minVal2   = min(min(min(NiiMap2_Masked)));

idx1      = find( NiiMap1_Masked > (MinRange*minVal1) & NiiMap1_Masked < (MaxRange*maxVal1) );
idx2      = find( NiiMap2_Masked > (MinRange*minVal2) & NiiMap2_Masked < (MaxRange*maxVal2) );
final_idx = intersect(idx1, idx2);

msk1_idx  = find(Mask1(:, :, SliceNum)~=0);
msk2_idx  = find(Mask2(:, :, SliceNum)~=0);
msk3_idx  = find(Mask3(:, :, SliceNum)~=0);
msk4_idx  = find(Mask4(:, :, SliceNum)~=0);

if Remove_Extreme_Vals
    not_max_idx_1 = find(NiiMap1_Masked ~= max(max(NiiMap1_Masked)));
    not_max_idx_2 = find(NiiMap2_Masked ~= max(max(NiiMap2_Masked)));
    not_min_idx_1 = find(NiiMap1_Masked ~= min(min(NiiMap1_Masked)));
    not_min_idx_2 = find(NiiMap2_Masked ~= min(min(NiiMap2_Masked)));
    
    new_idx       = intersect(not_max_idx_1, not_max_idx_2);
    new_idx       = intersect(new_idx, not_min_idx_1);
    new_idx       = intersect(new_idx, not_min_idx_2);
    final_idx     = intersect(new_idx,final_idx);
else
    new_idx       = final_idx;
end

vector1   = NiiMap1_Masked( final_idx);
vector2   = NiiMap2_Masked( final_idx);

%[R, PValue] = corrplot([vector1  vector2], 'testR','on');
%title('Total corr. plot');
%[Rval, Pval] = corrplot([vector1  vector2], 'testR','on');

display(['-I- Displaying ' num2str(length(vector1)) ' voxels.']);

idx1 = intersect(new_idx,msk1_idx);
idx1 = intersect(idx1,final_idx);
idx2 = intersect(new_idx,msk2_idx);
idx2 = intersect(idx2,final_idx);
idx3 = intersect(new_idx,msk3_idx);
idx3 = intersect(idx3,final_idx);
idx4 = intersect(new_idx,msk4_idx);
idx4 = intersect(idx4,final_idx);

% All together
if Remove_BV
    idx5 = unique([idx1 ;idx2 ;idx3]);
else
    idx5 = unique([idx1 ;idx2 ;idx3 ;idx4]);
end

Mask_Names{5} = {'All'};

title([char(Mask_Names{1}) ' corr. plot']);

Pval = zeros(5,1);
Rval = zeros(5,1);

for i = 1 : 5
    eval(['[Rmat, Pmat] = corrplot([NiiMap1_Masked( idx' num2str(i) ' )  NiiMap2_Masked( idx' num2str(i) ' )], ''testR'',''on'');']);
    
    if Close_CorrPlot_Fig
        h = gcf;
        close(h);
    end
    Rval(i) = Rmat(1,2);
    Pval(i) = Pmat(1,2);
    
end

fig_num = figure;
h1 = scatter(NiiMap1_Masked( idx1 ), NiiMap2_Masked( idx1 ), 'r');
hold on;
h2 = scatter(NiiMap1_Masked( idx2 ), NiiMap2_Masked( idx2 ), 'g');
h3 = scatter(NiiMap1_Masked( idx3 ), NiiMap2_Masked( idx3 ), 'b');
if ~Remove_BV
    h4 = scatter(NiiMap1_Masked( idx4 ), NiiMap2_Masked( idx4 ), 'c');
end
hold off;

if Remove_BV
    
    
    for i = 1 : 3
        eval(['legened' num2str(i) '_name = '' ' char(Mask_Names{i}) ' P= ' num2str(Pval(i),'%.2f') ' R= ' num2str(Rval(i),'%.2f') ''' ;']);
    end
    
    [hleg, hobj] = legend( [h1 h2 h3], legened1_name, legened2_name,legened3_name, 'Location','NorthWest');
    
else
    
    for i = 1 : 4
        eval(['legened' num2str(i) '_name = '' ' char(Mask_Names{i}) ' P= ' num2str(Pval(i),'%.2f') ' R= ' num2str(Rval(i),'%.2f') ''' ;']);
    end
    
    [hleg, hobj] = legend( [h1 h2 h3 h4 ], legened1_name, legened2_name,legened3_name, legened4_name, 'Location','NorthWest');
    
    
end






textobj      = findobj(hobj, 'type', 'text');
set(textobj, 'Interpreter', 'latex', 'fontsize', font_legend);

%title_string = sprintf(['Error Stats - ' param_name '. Best: ' num2str(best_error, '%.2f') '+-' num2str(best_std_or_sem, '%.2f') ]);
title_string = sprintf(['DCE-DSC Corr. P_{total}: ' num2str(Pval(end),'%.2f') ' R_{total}: ' num2str(Rval(end),'%.2f')]);

title(title_string,'FontSize',font_size,'FontWeight','bold');
xlabel(['DCE Norm. ' map_name ' values'],'FontSize',font_size,'FontWeight','bold');
ylabel(['DSC Norm. ' map_name ' values'],'FontSize',font_size,'FontWeight','bold');

saveas(fig_num, [OutPath filesep Subj_Name '_' map_name '_Correlation_Comparison_Values.jpg'])

%% Figure to put in article
if Scatter_no_legend
    
    fig_num = figure;
    h1 = scatter(NiiMap1_Masked( idx1 ), NiiMap2_Masked( idx1 ), 'or');
    hold on;
    h2 = scatter(NiiMap1_Masked( idx2 ), NiiMap2_Masked( idx2 ), '+g');
    h3 = scatter(NiiMap1_Masked( idx3 ), NiiMap2_Masked( idx3 ), '^b');
    hold off;
    
    for i = 1 : 3
        eval(['legened' num2str(i) '_name = '' ' char(Mask_Names{i}) ''' ;']);
    end
    [hleg, hobj] = legend( [h1 h2 h3], legened1_name, legened2_name,legened3_name, 'Location','NorthWest');
    textobj      = findobj(hobj, 'type', 'text');
    set(textobj, 'Interpreter', 'latex', 'fontsize', font_legend);
    
    title_string = sprintf(['DCE-DSC Scatter Plot']);
    
    title(title_string,'FontSize',font_size,'FontWeight','bold');
    xlabel(['DCE Norm. ' map_name ' Values'],'FontSize',font_size,'FontWeight','bold');
    ylabel(['DSC Norm. ' map_name ' Values'],'FontSize',font_size,'FontWeight','bold');
    
    set(gca,'FontSize',font_size,'FontWeight','bold');
    
    % Add Least Squares Line
    %lsline;
    
    saveas(fig_num, [OutPath filesep Subj_Name '_' map_name '_Correlation_Comparison.jpg'])
    
end

end



