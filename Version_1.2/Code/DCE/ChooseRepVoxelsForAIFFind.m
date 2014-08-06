function [CVI BinCVI Bin2CVI]=ChooseRepVoxelsForAIFFind(Msk3D,CTC2DOrig,BolusStart,nSets1,nSets2,nPerSet)

Msk2=Msk3D;

% Erodes each of the slices
for i=1:size(Msk3D,3)
    Msk2(:,:,i)=imerode(squeeze(Msk2(:,:,i)),strel('disk',4));
end

% Get mask's indices of original 3d mask
% Msk1Idxs=find(Msk3D);
% Get mask's indices of eroded 3d mask (with good indices of the original mask)
% GoodIdxsM=find(Msk2(Msk1Idxs));

% Out of the eroded mask, keep indices of the original 3D mask
GoodIdxsM=find(Msk2(Msk3D));

% Devide each slice by the maximum
rCTC2D=CTC2DOrig ./ repmat( max(CTC2DOrig,[],2) , [1 size(CTC2DOrig,2) ] );

% Mask  slices with negative ctc value
GoodIdxsR=find( sum( rCTC2D(:,BolusStart+2:end) < -0.001,2 ) < 3 );

% Intersect the eroded mask and the mask that checked for negative values
GoodIdxs=intersect(GoodIdxsM,GoodIdxsR);

% Get the masked CTC2D
CTC2D=CTC2DOrig(GoodIdxs,:);

%% Estimate noise for each voxel according to median C(t) divided by max C(t)

rmadCTC2D=EstimateNoise(CTC2D);
% adCTC2D=abs(diff(CTC2D,[],2));
% madCTC2D=median(adCTC2D,2);
% rmadCTC2D=madCTC2D./max(CTC2D,[],2);
% BadNoise=find(rmadCTC2D<0);

%% Auto points
MaxAroundBolus=max(CTC2D(:,(BolusStart-1):(BolusStart+3)),[],2);
MaxAtEnd=max(CTC2D(:,(end-5):end),[],2);

% ASK GILAD - Seems like he is doing here the same calculation as in DCET1_PKf.m
% ANSWER - One of them is redundant. Fixed in newer version.

%%
Stuff1=MaxAroundBolus;
Stuff2=MaxAroundBolus./MaxAtEnd;

S=sort(Stuff1);
N=size(CTC2D,1);
Vals=linspace(1,N,nSets1+1);
Edges=[S(floor(Vals(1:end-1))); max(S)+1];
[H,Bin] = histc(Stuff1,Edges);

% S=sort(MaxAtEnd);
% S=sort(MaxAtEnd./MaxAroundBolus);
% N=size(CTC2D,1);
% Vals=linspace(1,N,nSets2);
% Edges=S(floor(Vals));
% [H,Bin2] = histc(MaxAtEnd,Edges);
Bin2=Bin*0;
for i=1:nSets1
    VV=Stuff2(Bin==i);
    S=sort(VV);
    N=numel(S);
    Vals=linspace(1,N,nSets2+1);
    Edges=[S(floor(Vals(1:end-1))); max(S)+1];
    [H,Bin2C] = histc(VV,Edges);
    Bin2(Bin==i)=Bin2C;
end

Clust=(Bin-1)*nSets2+Bin2;
Reps=[];
Noises=[];

for i=1:max(Clust)
    CurI=find(Clust==i);
    Ns(i)=numel(CurI);
%     if(numel(CurI)<100)
%         continue;
%     end
%     [S Ord]=sort(rmadCTC2D(CurI));
    Ord=1:numel(S);
    Reps=[Reps CurI(Ord(1:nPerSet)')'];
    Noises=[Noises S(1:nPerSet)'];
end
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

% nChosen=numel(CVI);
% figure(303);clf;
% for i=1:nChosen
%     gsubplot(nChosen,i);
%     plot(CTC2DOrig(CVI(i),:),'*');
% end

a=5;