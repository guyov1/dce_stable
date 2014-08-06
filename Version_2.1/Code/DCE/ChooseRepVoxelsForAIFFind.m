function [CVI BinCVI Bin2CVI DataToFit]=ChooseRepVoxelsForAIFFind(Msk3D,MskE,CTC2DOrig,BolusStart,nSets1,nSets2,nPerSet)

% Msk2=Msk3D;
% MskE=Msk2;
% FBrainMask=bwfillHoles3Dby2D(Msk2);

% % Erodes each of the slices
% for i=1:size(Msk3D,3)
%     %figure;
%     %subplot(2,1,1);imshow(Msk2(:,:,i));title('before');
% %     Msk2(:,:,i)=imerode(squeeze(Msk2(:,:,i)),strel('disk',10));
%     MskE(:,:,i)=imerode(squeeze(FBrainMask(:,:,i)),strel('disk',16));
%     %subplot(2,1,2);imshow(Msk2(:,:,i));title('After');
% end
Msk2=Msk3D & MskE;

nVoxels=size(CTC2DOrig,1);

% Get mask's indices of original 3d mask
% Msk1Idxs=find(Msk3D);
% Get mask's indices of eroded 3d mask (with good indices of the original mask)
% GoodIdxsM=find(Msk2(Msk1Idxs));

% Out of the eroded mask, keep indices of the original 3D mask
GoodIdxsM=find(Msk2(Msk3D));

% Devide each slice by the maximum
rCTC2D=CTC2DOrig ./ repmat( abs(max(CTC2DOrig,[],2)) , [1 size(CTC2DOrig,2) ] );

% Mask  slices with negative ctc value
% GoodIdxsR=find( sum( rCTC2D(:,BolusStart+2:end) < -0.001,2 ) < 3 );
GoodIdxsRa=find( sum( rCTC2D(:,max(1,BolusStart-5):end) < -0.001,2 ) < 3);
GoodIdxsRb=find( sum( rCTC2D(:,max(1,BolusStart-5):end) < -0.05,2 ) < 1);
if(numel(GoodIdxsRb)<5)
    GoodIdxsRb=1:nVoxels;
end
% GoodIdxsR=find( sum( rCTC2D(:,max(1,BolusStart-5):end) < -0.001,2 ) < 3 & sum( rCTC2D(:,max(1,BolusStart-5):end) < -0.05,2 ) < 1);
GoodIdxsR=intersect(GoodIdxsRa,GoodIdxsRb);

% Remove the smallet ones
TmpA=max(CTC2DOrig,[],2);
Tmp=sort(TmpA);
Thresh=Tmp(max(1,floor(numel(Tmp)*0.0)));
GoodIdxsS=find(TmpA>Thresh);

% Remove the bigget ones
Tmp=sort(Tmp,'descend');
Thresh=median(Tmp(1:min(numel(Tmp),1000)))*20;
GoodIdxsB=find(TmpA<Thresh);

% Intersect the eroded mask and the mask that checked for negative values
GoodIdxs=intersect(intersect(GoodIdxsM,GoodIdxsR),intersect(GoodIdxsS,GoodIdxsB));

% Get the masked CTC2D
CTC2D=CTC2DOrig(GoodIdxs,:);

if(isempty(GoodIdxs))
    CVI=[];
    return;
end
%% Auto points
MaxAroundBolus = max(CTC2D(:,(BolusStart-1):(BolusStart+3)),[],2);
MaxAtEnd       = max(CTC2D(:,(end-5):end),[],2);

% ASK GILAD - Seems like he is doing here the same calculation as in DCET1_PKf.m
% ANSWER - One of them is redundant. Fixed in newer version.

%% Get temporal interesting volumes around the bolus time and end time (it will characterize each C(t)).
% Each voxel has a c(t). The interesting characteristics of each voxel 
% are the c(t) at bolus arrival and c(t) at the end of test (to know if we have enhancement).
% We used those values to cluster all the voxels so eventually we could work on some representatives.
Stuff1  = MaxAroundBolus;
Stuff2  = MaxAroundBolus./MaxAtEnd;

% Sort the volumes around bolus time according to maximal value
S       = sort(Stuff1);
% Number of masked pixels in each temproal volume
N       = size(CTC2D,1);
% Line space 5 points between 1 and the number of masked pixels in each volume
Vals    = linspace(1,N,nSets1+1);
% Get the 5 values of the sorted maximal value pixels
Edges   = [S(floor(Vals(1:end-1))); max(S)+1];
% Histogram count for all the pixels around the bolus according to maximal
% and minimal value (1) divided to 5 bins
[H,Bin] = histc(Stuff1,Edges);

% S=sort(MaxAtEnd);
% S=sort(MaxAtEnd./MaxAroundBolus);
% N=size(CTC2D,1);
% Vals=linspace(1,N,nSets2);
% Edges=S(floor(Vals));
% [H,Bin2] = histc(MaxAtEnd,Edges);
Bin2=Bin*0;
% Do the same sorting and binning for values close to the end of experiment
for i=1:nSets1
    VV          = Stuff2(Bin==i);
    if(isempty(VV))
        continue;
    end
    S           = sort(VV);
    N           = numel(S);
    Vals        = linspace(1,N,nSets2+1);
    if(numel(VV)==1)
        Edges       = [S(floor(Vals(1:end-1)))'; max(S)+1];
    else
        Edges       = [S(floor(Vals(1:end-1))); max(S)+1];
    end
    [H,Bin2C]   = histc(VV,Edges);
    Bin2(Bin==i)= Bin2C;
end

% Doing a global clustering (both for bolus time and end time). 
% If one point has index (2,3), it is actually 2*5 + 3.
Clust=(Bin-1)*nSets2+Bin2;
Reps   = [];
Noises = [];

rmadCTC2D=EstimateNoise(CTC2D);
%% Get representing voxels according to lowest noise

% For each cluster, get the least noisiest 2 voxels (and their index)
for i=1:max(Clust)
    % Indices of all voxels in cluster i
    CurI  = find(Clust==i);
    % Number of voxels in that cluster
    Ns(i) = numel(CurI);
    
    % If there are not enough voxels in this cluster (<100), continue
    if(size(CTC2D,1)>50 && numel(CurI)<10)
        continue;
    end
    if(isempty(CurI))
        continue;
    end
    
    % Get a sorted list (S) of the noise medians of all the voxels from the
    % specific cluster. Get the original indexing (Ord) as well.
    [S Ord]= sort(rmadCTC2D(CurI));
    
%     Ord=1:numel(S);

    % Get the voxels indices with the loweset noise values (2 smallest)
    Reps   = [Reps CurI(Ord(1:min(numel(S),nPerSet))')'];
    % Get the noise values
    Noises = [Noises S(1:min(numel(S),nPerSet))'];    
end

% Sort the noises from all clusters (+ original indices)
% [S Ord]=sort(Noises);
Ord=1:numel(Noises);
% Rep2D=CTC2D(Reps(Ord),:)./repmat(max(CTC2D(Reps(Ord),:),[],2),[1 size(CTC2D,2)]);
% Good=find(sum(Rep2D(:,BolusStart+2:end)<-0.001,2)<3);

% Rep2D=Rep2D(Good,:);
% CVI=Reps(Ord(Good));
CVI=Reps(Ord);
%%
BinCVI=Bin(CVI);
Bin2CVI=Bin2(CVI);

CVI=GoodIdxs(CVI);
% CVI=setdiff(CVI,BadNoise);
DataToFit=CTC2DOrig(CVI,:);
% nChosen=numel(CVI);
% figure(303);clf;
% for i=1:nChosen
%     gsubplot(nChosen,i);
%     plot(CTC2DOrig(CVI(i),:),'*');
% end

a=5;