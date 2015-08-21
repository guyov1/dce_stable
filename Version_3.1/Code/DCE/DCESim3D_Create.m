BaseP='\\fmri-t9\users\Moran\DCE\ForSim\';
D=dir(BaseP);
D=D([D.isdir]);
D=D(3:end);
D=D(~strhas({D.name},'SmVl'));
D=D(~strhas({D.name},'LiHa'));
%%
for p=1:numel(D)
    try
        disp('a');
        WorkingP=[BaseP D(p).name filesep];
        USStr='';
        SimFN=[WorkingP 'Sim.mat'];
        if(exist(SimFN,'file'))
            disp([SimFN ' already exists, skipping.']);
            continue;
        end
        SimPKMFN=[WorkingP 'SimPKM' USStr '.mat'];
        SimPKM3DFN=[WorkingP 'SimPKM3D' USStr '.mat'];
        % Load CTCs and get noise per pixel
        CTCFN=[WorkingP 'AfterCTC' '.mat'];
        a=load(CTCFN);
        NoisePerPixel=median(abs(diff(a.CTC2DBigGood,[],2)),2);
        % Load AIF
        disp('Load AIF');
        PKMFN=[WorkingP 'PKM' USStr '.mat'];
        b=load(PKMFN,'OutAIFParam');
        PrepareFN=[WorkingP 'AfterPrepare4D.mat'];
        c=load(PrepareFN,'TimeBetweenDCEVols','BolusStart');
        nSVols=size(a.CTC2D,2);
        TimeBetweenDCEVolsMin=c.TimeBetweenDCEVols/60;
        InterpolationFactor=ceil(c.TimeBetweenDCEVols);
        SampleTs=((1:nSVols)-1)*TimeBetweenDCEVolsMin;
        
        HInterpolationFactor=ceil(InterpolationFactor*2);
        Hdt=TimeBetweenDCEVolsMin/HInterpolationFactor;
        HSampleTs=0:Hdt:SampleTs(end);
        
        ThreeSec=ceil(3/(Hdt*60))*2;
        TDif=-Hdt*ThreeSec:Hdt:Hdt*ThreeSec;
        nTDif=numel(TDif);
        
        AIF_Parker8t=@(x,t) AIF_Parkerg2( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8))*x(2);
        HAIFOrig=AIF_Parker8t(b.OutAIFParam,HSampleTs);
        % Load PKs
        disp('Load PKs');
        PKM3DFN=[WorkingP 'PKM3D' USStr '.mat'];
        d=load(PKM3DFN,'PKs');
        % Distort BAT and PKs
        N=size(a.CTC2DBigGood,1);
        tmpBATIs=d.PKs(:,1);
        BNans=isnan(tmpBATIs);
        tmpBATIs(BNans)=1;
        SBATs=TDif(tmpBATIs)'+rand(N,1)/60-(1/120);
        SBATs(BNans)=rand(sum(BNans),1)/60-(1/120);
%    1        2         3        4          5     6        7    8   9    10
% BATfinal VpFinal KtransFinal Kepfinal VeFinal RSSFinal RSS0 RSS1 RSS2 RSS3
        
        SPKs=max(d.PKs(:,[1 4 2 3]).*(randn(size(d.PKs(:,[1 4 2 3])))*0.1+1),0);
%       
MinVp=0;
MaxVp=1;
MinVe=0.05;
MaxVe=10; % Account for AIF maxAmp too big -> ktrans small
MaxKtrans=2.2;
MaxKep=MaxKtrans/MinVe;
% BAT Kep Vp Ktrans Ve RMS fVal f0Val


SPKs(:,2)=min(MaxKep,SPKs(:,2));
SPKs(:,3)=min(MaxVp,SPKs(:,3));
SPKs(:,4)=min(MaxKtrans,SPKs(:,4));

        MinFirstBolusSig=2; % seconds
        LB=[0.1 0     MinFirstBolusSig/60 0.1     0   0.1  0.1 0]';
        UB=[10  1.5   0.25                10      1   3    2   0.3]';
% 
        NewParams=(b.OutAIFParam).*(randn(size(b.OutAIFParam))*0.1+1);
%         NewParams=min(UB',max(LB',NewParams));
        HAIF=AIF_Parker8t(NewParams,HSampleTs);
%         %% Create Sims
%         disp('Create Sims');
%         disp(WorkingP);
        HConvIdxM=CreateConvIdxMFromSampleTs(numel(HSampleTs));
        HTriB=HConvIdxM>0;
        HConvIdxMTriB=HConvIdxM(HTriB);
%         Sims=a.CTC2DBigGood*0;
%         
        NAtATime=10000;
        Sims=nan(N,numel(SampleTs));
        disp(['There are ' num2str(N) ' voxels to simulate ' WorkingP]);
        before=now;
        for i=1:NAtATime:N
            tic
            IdxsA=i:min(N,i+NAtATime-1);
            AllZeros=find(sum(SPKs(IdxsA,3:4).^2,2)==0);
            Idxs=setdiff(IdxsA,AllZeros);
            HConvd2=DCECostFuncgrT1ForConv(HAIF',SPKs(Idxs,2),HSampleTs,HConvIdxMTriB,HTriB);
            for j=1:numel(Idxs)
                SAIF=interp1(HSampleTs,HAIF,SampleTs+SBATs(Idxs(j)),'linear','extrap');
                SHConvd2=interp1(HSampleTs,HConvd2(j,:)',SampleTs+SBATs(Idxs(j)),'linear','extrap')';
                Regressors=[SAIF; SHConvd2'];
                Sims(Idxs(j),:)=((Regressors')*(SPKs(Idxs(j),[3 4])'));
            end
            Sims(AllZeros,:)=0;
            t=toc;
            TimeFromStart=now-before;
            WillFinishAt=before+TimeFromStart*N/Idxs(end);
            disp(['Simulating ' num2str(Idxs(1)) ':' num2str(Idxs(end)) ' took ' num2str(t) 's. Will finish around ' datestr(WillFinishAt)]);
        end
        disp(['Finished sim ' WorkingP]);
        % BAT Kep Vp Ktrans Ve RMS fVal f0Val
        % figure;plot(SampleTs,a.CTC2DBigGood(Idxs(9),:),'k*',SampleTs,Sims(Idxs(9),:),'b')
        % Add noise
        NSims=Sims+randn(size(Sims)).*repmat(NoisePerPixel,[1 nSVols]);
        % i=200;figure;plot(SampleTs,NSims(i,:),'*',SampleTs,Sims(i,:),'k-',SampleTs,a.CTC2DBigGood(i,:),'r*')
        save(SimFN,'SBATs','SPKs','Sims','NSims','NewParams','HAIF','HAIFOrig');
    catch ME
        MEs{p}=ME;
        disp(['Error in ' num2str(p) ' ' D(p).name]);
    end
end