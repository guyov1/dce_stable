BaseP='\\fmri-t9\users\Moran\DCE\ForSim\';
D=dir(BaseP);
D=D([D.isdir]);
D=D(3:end);
D=D(~strhas({D.name},'SmVl'));
D=D(~strhas({D.name},'LiHa'));
%%
for p=1:numel(D)
    try
        UseBAT=false;
        if(UseBAT)
            BATStr='WithBAT';
        else
            BATStr='WithoutBAT';
        end
        disp('a');
        WorkingP=[BaseP D(p).name filesep];
        USStr='';
        SimFN=[WorkingP 'Sim.mat'];
        SimPKMFN=[WorkingP 'SimPKM' USStr '.mat'];
        SimPKM3DFN=[WorkingP 'SimPKM3D' USStr '_' BATStr '.mat'];
        if(exist(SimPKM3DFN,'file'))
            disp([SimPKM3DFN ' already exists, skipping.']);
            continue;
        end
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
        if(~UseBAT)
            ThreeSec=0;
        end
        
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
%         SBATs=TDif(d.PKs(:,1))'+rand(N,1)/60-(1/120);
%         SPKs=max(d.PKs.*(randn(size(d.PKs))*0.1+1),0);
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
%         NewParams=(b.OutAIFParam).*(randn(size(b.OutAIFParam))*0.1+1);
%         NewParams=min(UB',max(LB',NewParams));
%         HAIF=AIF_Parker8t(NewParams,HSampleTs);
%         %% Create Sims
%         disp('Create Sims');
%         disp(WorkingP);
        HConvIdxM=CreateConvIdxMFromSampleTs(numel(HSampleTs));
        HTriB=HConvIdxM>0;
        HConvIdxMTriB=HConvIdxM(HTriB);
%         Sims=a.CTC2DBigGood*0;
%         
        NAtATime=10000;
        load(SimFN,'SBATs','SPKs','Sims','NSims','NewParams','HAIF','HAIFOrig');
        % Calc AIF
        
%         F1=a.F1(a.MskCTCGood);
%         B3=ismember(F1,a.F2);
% 
%         CTC2DA=NSims(B3(:),:);
%         
%         ImagB=any(imag(CTC2DA)~=0,2) | any(isnan(CTC2DA),2);
%         CTC2D=CTC2DA(~ImagB,:);
%         
%         Msk2=a.DBrainMask;
%         Msk2(Msk2)=a.MskCTCGood;
%         Msk2(Msk2)=B3;
%         Msk2(Msk2)=~ImagB;
%         
%         Idx3D=NaN(size(Msk2));
%         Idx3D(Msk2)=1:sumn(Msk2);
%         [CVI BinCVI Bin2CVI]=ChooseRepVoxelsForAIFFind(Idx3D>0,CTC2D,a.BolusStart,10,10,1);
%         MaxAmp=getMaxAmp(CTC2D);
%         
%         Options=struct('SubSecRes',{[]},'MaxTDif',{[]});
%         
%         adCTC2D=abs(diff(CTC2D,[],2));
%         madCTC2D=median(adCTC2D,2);
%         rmadCTC2D=madCTC2D./max(CTC2D,[],2);
%         
%         DataNoise=rmadCTC2D(CVI);
%         DataToFit=CTC2D(CVI,:);
%         
%         
%         T1=0.17046;A1=0.809;sig1=0.0563;sig2=0.132;A2=0.330;T2=0.365;alpha=1.050;beta=0.1685;s=38.078;tau=0.483;
%         tauDelta=tau-T1;T2Delta=T2-T1;
%         % AIFParamsA=[T1 A1 sig1];
%         T1=1;
%         
%         ParamAIFCoeff=[0.8,25,3;... % First bolus time
%             0.25,4,1;... % First bolus magnitude
%             0.25,10,2;... % First bolus width
%             0.5,10,1;... % 2nd bolus time
%             0.25,4,1;... % 2nd bolus magnitude
%             0.25,4,1;... % 2nd bolus width
%             0.5,2,1;... % alpha
%             0.5,2,1;... % beta
%             0.5,2,1;... % s
%             0.5,5,1]; % tauDelta
%         
%         ParamAIFCoeff(1,3)=a.BolusStart*a.TimeBetweenDCEVols/60;
%         ParamAIFCoeff(3,:)=[1 5 3];
%         ParamAIFCoeff(1,1)=0.1;
%         ParamAIFCoeff(1,2)=SampleTs(end)/2;
%         ParamAIFCoeff(5,1)=0;
%         ParamAIFCoeff(7,:)=[0.2 2 0.4];
%         ParamAIFCoeff(8,:)=[0 0.3 beta];
%         % ParamAIFCoeff(7,2)=20;
%         % ParamAIFCoeff(8,2)=20;
%         ParamAIFCoeff(10,1)=0;
%         ParamAIFCoeff(10,3)=0.5;
%         ParamAIFCoeff(2,2)=10;
%         ParamAIFCoeff(6,:)=[1 5 3];
%         
%         nAIFParams=size(ParamAIFCoeff,1);
%         
%         SimAIFFinderFN=[WorkingP 'SimAIFFindData' USStr '.mat'];
%         
%         [PKOut OutAIFParam]=AIFTryf(WorkingP,DataToFit,MaxAmp,DataNoise,ParamAIFCoeff,nSVols,a.TimeBetweenDCEVols,Options,false,SimAIFFinderFN);
%         save(SimPKMFN,'OutAIFParam','DataToFit');
%         
%         EAIF=AIF_Parker8t(OutAIFParam,HSampleTs);
        
        EAIF=HAIF;
%         figure(2000);clf;plot(HSampleTs,HAIF,'k','LineWidth',4);hold on;plot(HSampleTs,EAIF,'m');
        % Calc PKs
        CHAIF=cumtrapz(HSampleTs,EAIF);
        
        SAIF=zeros([nTDif numel(SampleTs)]);
        CSAIF=zeros([nTDif numel(SampleTs)]);
        for i=1:nTDif
            SAIF(i,:)=interp1(HSampleTs,EAIF,SampleTs+TDif(i),'linear','extrap');
            CSAIF(i,:)=interp1(HSampleTs,CHAIF,SampleTs+TDif(i),'linear','extrap');
        end
        
        MskCTCGood3D=a.DBrainMask;
        MskCTCGood3D(MskCTCGood3D)=a.MskCTCGood;
        %
        PKs=NaN(N,28);
        NAtATime=5000;
        disp(['There are ' num2str(N) ' voxels to compute ' WorkingP]);
        before=now;
        for i=1:NAtATime:N
            tic
            CurIs=i:min(N,i+NAtATime-1);
%             PKs(CurIs,:) = FindPKBATgAIFMuraseF_TProb(NSims(CurIs,:),SAIF,SampleTs,CSAIF);
                                                            %DataToFit,SAIF,SampleTs,CumSumAIFdTSamp,WPerCTC)
            PKs(CurIs,:) = FindPKBATgAIFMuraseF4Models_TProb(NSims(CurIs,:),SAIF,SampleTs,CSAIF);
            t=toc;
            TimeFromStart=now-before;
            WillFinishAt=before+TimeFromStart*N/CurIs(end);
            disp(['Calculating ' num2str(CurIs(1)) ':' num2str(CurIs(end)) ' took ' num2str(t) 's. Will finish around ' datestr(WillFinishAt)]);
        end
        disp(['Finished sim PK estimation ' WorkingP]);
        save(SimPKM3DFN,'PKs');
        %
%         EBATs=TDif(PKs(:,1));
    catch ME
        MEs{p}=ME;
        disp(['Error analyze in ' num2str(p) ' ' BATStr ' ' D(p).name]);
    end
end