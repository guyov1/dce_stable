function DCET1_RelaxForSubjectf(DCEInfos,WorkingP,CalcForce,WhichMean,Options,MainInfo)

AddToLog(WorkingP,'b_00','\\subsection*{Relaxometry}');

% DoN3,DoGlobal,DoDeviations
RelaxFN=[WorkingP 'Relax.mat'];
%% Create the Relaxometry directory and needed initial structs in case it does not exist already
DCERelaxP=[WorkingP 'Relaxometry' filesep];

% Create relaxometry directory
mkdir(DCERelaxP);


%% Get the number of flip angle sets
%   NOTE - I assume all the sets have the same angle scans (repeated for each set)

% First angle used
first_angle = DCEInfos(1).FlipAngle;
% All indices with same angle
index_new_sets = find([DCEInfos.FlipAngle]==first_angle);
% Number of sets
num_FAs_sets = numel(index_new_sets);

% Create a copy of DCE Infos from which we pick every iteration the relevant infos
Copy_DCEInfos = DCEInfos;

%% Repeat the T1 process for every set of angles
for j =1:num_FAs_sets
    
    % Create a directory for each set
    Set_Dir =[DCERelaxP 'Set_Num_' sprintf('%d',j) filesep];
    mkdir(Set_Dir);
    % Set the working directory to be the same
    RWorkingP=Set_Dir;
    
    % Get the relevant DCE infos to work on
    if (j ~= num_FAs_sets)  % if not last iteration
        first_angle_to_work = index_new_sets(j);
        last_angle_to_work   = index_new_sets(j+1) -1;
        
        % If we have less than 2 angles to work with, we can not create T1
        % map -> Error!
        if (last_angle_to_work - first_angle_to_work + 1 < 2)
            helpdlg('Error in Data!!! One of the DESPOT T1 data has less than 2 flip angles. Can not create T1 map this way!');
            error('Error in Data!!! One of the DESPOT T1 data has less than 2 flip angles. Can not create T1 map this way!');
            
        end
        
        DCEInfos = Copy_DCEInfos(first_angle_to_work :last_angle_to_work);
    else % last iteration treatment
        first_angle_to_work = index_new_sets(j);
        
        % If we have less than 2 angles to work with, we can not create T1
        % map -> Error!
        if ( size(Copy_DCEInfos,2) - first_angle_to_work + 1 < 2)
            helpdlg('Error in Data!!! One of the DESPOT T1 data has less than 2 flip angles. Can not create T1 map this way!');
            error('Error in Data!!! One of the DESPOT T1 data has less than 2 flip angles. Can not create T1 map this way!');
            
        end
        
        DCEInfos = Copy_DCEInfos(first_angle_to_work :end);
    end
    
    % If calculated T1 maps already, continute to next iteration
    % If we got to last iteration, and T1 map was calculated, return
    % DoDeviations - with N3 cleaner
    %     if (DoDeviations)
    if Options.ExtractFAs
        if Options.UseN3OnT1
            T1_Map_Name = 'T13DNFA_N3k.nii';
        else
            T1_Map_Name = 'T13DNFA.nii';
        end
    else
        if Options.UseN3OnT1
            T1_Map_Name = 'T13DOFA_N3k.nii';
        else
            T1_Map_Name = 'T13DOFA.nii';
        end
    end
    
    if(exist([RWorkingP T1_Map_Name],'file') && ~CalcForce)
        
        disp('DCET1_RelaxForSubjectf Already computed');
        
        % If last iteration, return, else, continue loop
        if (j == num_FAs_sets)
            return;
        else
            continue;
        end
    end
    
    disp('Starting DCET1_RelaxForSubjectf');
    
    % Get repetetion time (TRs) of all the scans
    AllInfos=DCEInfos;
    if(j==1 && Options.IncludeMainInT1)
        AllInfos=[AllInfos MainInfo];
    end
    
    % Get flip angle (FAs) of all the scans
    FAsF=[AllInfos.FlipAngle];
    TRsF=[AllInfos.RepetitionTime];
    TEsF=[AllInfos.EchoTime];
    R1sF=[AllInfos.R1];
    R2sF=[AllInfos.R2];
    TGsF=[AllInfos.TG];
    
    % Get average Time Stamp of every scan
    Time_Stamps =   {DCEInfos.SeriesDateTime};
    size_of_list = size(Time_Stamps,2);
    % Get the middle index to ve the time stamp
    middle_index = ceil(size_of_list/2);
    Mid_Time_Stamp = [Time_Stamps(middle_index)];
    
    % Create a text file which will hold the average time of the scan
    file_name = [RWorkingP 'Middle_Time_Stamp.mat'];
    save (file_name,'Mid_Time_Stamp')
    
    
    %% Base
    CMeanFN=[WorkingP 'CopyDCEMean.nii'];
    % Get all FAs information
    CurOtherFAsInfos=DCEInfos;
    % Number of FAs infos
    nOtherFAs=numel(CurOtherFAsInfos);
    % Prepare struct list of Base/Co-Registered and Mat
    OtherFACrgFNs=cell(1,nOtherFAs);
    OtherFABaseFNs=cell(1,nOtherFAs);
    OtherFAMatFNs=cell(1,nOtherFAs);
    
    %% Dicom2Nifti
    
    % Creating nifti for all original FAs
    for i=1:nOtherFAs
        
        % Create file name for nifti
        OtherFABaseFNs{i}=[Set_Dir 'OrigFA_' sprintf('%02d',CurOtherFAsInfos(i).FlipAngle) '_' sprintf('%02d',CurOtherFAsInfos(i).SeriesNumber) '.nii'];
        
        % Create temp nifti out of dicom data
        Tmp=gDicom2Nifti(CurOtherFAsInfos(i).Path,OtherFABaseFNs{i});
        
        % Same nifti file name concatenated with "_01"
        MainSubFN=[Set_Dir 'OrigFA_' sprintf('%02d',CurOtherFAsInfos(i).FlipAngle) '_' sprintf('%02d',CurOtherFAsInfos(i).SeriesNumber) '_01.nii'];
        
        % Not clear what the purpose of the following...
        if(Tmp>1 && Tmp<10 && exist(MainSubFN,'file'))
            
            movefile(MainSubFN,OtherFABaseFNs{i});
            delete([Set_Dir 'OrigFA_' sprintf('%02d',CurOtherFAsInfos(i).FlipAngle) '_' sprintf('%02d',CurOtherFAsInfos(i).SeriesNumber) '_*.nii']);
            
        end
        
        if(~exist(OtherFABaseFNs{i},'file'))
            error('DCET1_RelaxForSubjectf Dcm2Nii NotCreated');
        end
        
    end
    
    %% Coreg dest
    
    MedFA=median(FAsF);
    MedFAI=find(FAsF==MedFA,1);
    
    if(strcmp(WhichMean,'Mean 4D'))
        MeanFN=[WorkingP 'DCEMean.nii'];
    else
        if(strcmp(WhichMean,'Median angle'))
            MeanFN=OtherFABaseFNs{MedFAI};
        else
            MeanFN=WhichMean;
        end
    end
    delete(CMeanFN);
    try
        copyfile(MeanFN,CMeanFN,'f');
    catch MErr
    end
    
    if(~exist(CMeanFN,'file'))
        error('DCET1_RelaxForSubjectf: Problem copying CMeanFN');
    end
    
    %% Coreg
    for i=1:nOtherFAs
        OtherFAMatFNs{i}=CoregEstimate(OtherFABaseFNs{i},CMeanFN,CalcForce);
        OtherFACrgFNs{i}=CoregWrite(OtherFABaseFNs{i},OtherFAMatFNs{i},CalcForce,Set_Dir,false,CMeanFN);
    end
    delete(CMeanFN);
    disp('Relaxometry coregistration finished');
    
    %% Full brain mask for relaxometry
    PrepareFN=[WorkingP 'AfterPrepare4D.mat'];
    
    %load(PrepareFN,'BrainMask','C2','RepSli','MeanVol');
    load(PrepareFN,'BrainMask','RepSli','MeanVol');
    
    % If the brain mask have black holes in it, fill them
    FBrainMask=bwfillHoles3Dby2D(BrainMask);
    A=load_untouch_nii(OtherFACrgFNs{1});
    A.img=int16(FBrainMask);
    FMaskFN=[Set_Dir 'FBrainMsk.nii'];
    save_untouch_nii(A,FMaskFN);
    disp('FBrainMsk finished');
    
    %% Segmentation for masks
    se=strel('disk',4,8);
    EBrainMask=imerode(FBrainMask,se);
    
    %%
    AddToLog(WorkingP,'b_2d',['TRs:        ' num2str(TRsF,'% 8.2f')]);
    AddToLog(WorkingP,'b_2e',['TEs:        ' num2str(TEsF,'% 8.2f')]);
    AddToLog(WorkingP,'b_2f',['FAs:        ' num2str(FAsF,'% 8.2d')]);
    AddToLog(WorkingP,'b_2f2',['R1s:        ' num2str(R1sF,'% 8.2d')]);
    AddToLog(WorkingP,'b_2f3',['R2s:        ' num2str(R2sF,'% 8.2d')]);
    AddToLog(WorkingP,'b_2f4',['TGs:        ' num2str(TGsF,'% 8.2d')]);
    
    T2SApproximation=16;
    AddToLog(WorkingP,'b_2g',['Using $T_2^*$ approximation of ' num2str(T2SApproximation) 'ms']);
    T2SFactors=exp(-TEsF/T2SApproximation);
    AddToLog(WorkingP,'b_2h',['$T_2^*$ Factors: ' num2str(T2SFactors,'% 8.2f')]);
    T2SSigPerDif=100*(sqrt(max(T2SFactors)/min(T2SFactors))-1);
    AddToLog(WorkingP,'b_2i',['Estimated signal change because of echo time differences ' num2str(T2SSigPerDif) '\\%']);
    %%
    TRT1difF=100*(sqrt(max(TRsF)/min(TRsF))-1);
    AddToLog(WorkingP,'b_2j',['Estimated $T_1$ inaccuracy because of repetition time differences on all ' num2str(TRT1difF,'%2.2f') '\\%']);
    %%
    U=[-Inf unique(TRsF) Inf];
    MinGScr=Inf;
    for ii=1:numel(U)
        for jj=(ii+1):numel(U)
            CurGroupB=(TRsF>U(ii) & TRsF<U(jj));
            CurGroupF=find(CurGroupB);
            if(isempty(CurGroupF))
                continue;
            end
            CurSz=numel(CurGroupF);
            CurTRT1difF=(sqrt(max(TRsF(CurGroupF))/min(TRsF(CurGroupF)))-1);
            CurScr=1/CurSz+CurTRT1difF;
            if(CurScr<MinGScr)
                MinGScr=CurScr;
                BestGroup=CurGroupF;
            end
        end
    end
    AddToLog(WorkingP,'b_2l',['Chosen group to for FA extraction: ' num2str(BestGroup,'% d')]);
    %% Partition into groups, run on every group and compare results
    IncludeBaseline=Options.IncludeMainInT1;
    %     if(IncludeBaseline)
    %         if((exp(abs(log(TEsF(1)/TEsF(end))))-1)>0.2)
    %             IncludeBaseline=false;
    %             AddToLog(WorkingP,'b_2TE','TE difference, removing baseline from calculation');
    %         end
    %         if((exp(abs(log(TRsF(1)/TRsF(end))))-1)>0.3)
    %             IncludeBaseline=false;
    %             AddToLog(WorkingP,'b_2TR','TR difference, removing baseline from calculation');
    %         end
    %     end
    if(j==1 && IncludeBaseline)
        nOtherFAs=nOtherFAs+1;
        OtherFACrgFNs{nOtherFAs}=[WorkingP 'Baseline.nii'];
    end
    
    clear TRGroupC
    TRGroupC{1}=1:nOtherFAs;
    CurCGroup=1;
    
    clear FA4DF
    for i=1:nOtherFAs
        FA4DF(:,:,:,i)=loadniidata(OtherFACrgFNs{i});
    end
    WarningStatus=warning('off','MATLAB:rankDeficientMatrix');
    
    nFAGroups=numel(TRGroupC);
    ResNFAs=cell(nFAGroups,2);
    T1Res=cell(nFAGroups,3,2);
    
    TRGrpName=['[' regexprep(num2str(TRGroupC{CurCGroup}),'\D*','_') ']'];
    % RWorkingP=[BaseWorkingP TRGrpName filesep];
    MatFN=[RWorkingP 'NFARes.mat'];
    if(exist(MatFN,'file'))
        disp('Skipping');
        load(MatFN);
        %         continue;
    end
    
    TRGroup=TRGroupC{CurCGroup};
    Needed=TRGroup; % [1:4];
    FAs=FAsF(Needed);
    TRs=TRsF(Needed);
    TEs=TEsF(Needed);
    FA4D=FA4DF(:,:,:,Needed);
    nFAs=numel(FAs);
    SFA=size(FA4D);
    ANFAs=zeros(SFA(3),nFAs);
    ValMsk=all(FA4D>100,4);
    % FBrainMaskS=FBrainMask & ValMsk;
    
    %     CalcRelaxForVolNFA;
    if(~exist('ONFAs','var'))
        ONFAs=FAs;
    end
    %     load(MatFN,'ONFAs','FFA2DData','FBrainMaskS');
    %     load(MatFN,'ONFAs');
    
    FBrainMaskS=BrainMask;
    
    FFA2DData=Reshape4d22d(FA4D,FBrainMaskS>-1);
    FFA2DData=repMulti(FFA2DData,1./T2SFactors);
    %% Gain factors
    GainFactors=ones(1,nFAs);
    FFBrainMaskSF=find(FBrainMaskS);
    Tmp=FFA2DData(FFBrainMaskSF,:);
    Tmp=Tmp(getKrandomSamples(size(Tmp,1),1000),:);
    
    for tt=1:4
        for ii=1:nFAs
            Others=setdiff(1:nFAs,ii);
            %         Others=1:4;
            [RelaxCG{1:3}]=CalcT1byFAfw2(Tmp(:,Others)',FAs(Others),TRs(Others));
            SimByNomFA=SPGRfM(RelaxCG{1}',RelaxCG{2}',FAs(ii),TRs(ii));
            EstFac(ii)=median(Tmp(:,ii)./SimByNomFA);
        end
        AllEstGainFactors(tt,:)=EstFac;
        CurGainFactors=2.^round(log2(EstFac));
        AllGainFactors(tt,:)=CurGainFactors;
        Tmp=repMulti(Tmp,1./CurGainFactors);
    end
    NoUseCurrently_GainFactors=prod(AllGainFactors,1);
    Gains3=NoUseCurrently_GainFactors.*AllEstGainFactors(end,:);
    AddToLog(WorkingP,'b_3n1',['Estimated gains, not used but check:        ' num2str(AllEstGainFactors(1,:),'% 8.2f')]);
    AddToLog(WorkingP,'b_3n2',['Rounded gains, not used but check:        ' num2str(NoUseCurrently_GainFactors,'% 8.2f')],[],2);
    %     FFA2DData=repMulti(FFA2DData,1./GainFactors);
%     FFA2DData=repMulti(FFA2DData,1./Gains3);
    %%
    %         disp(['Baseline factor found: ' num2str(BestBaseFac)]);
    %         AddToLog(WorkingP,'b_Fac',['Baseline factor found: ' num2str(BestBaseFac)]);
    %         if(BestBaseFac<0.7 && BestBaseFac>0.3)
    %             disp('Changing to 0.5');
    %             Factors(end)=0.5;
    %         end
    %         FA4D(:,:,:,end)=FA4D(:,:,:,end)/Factors(end);
    %     end
    %% Calculate T1 maps out of the FAs registered images
    % The next line can take 1-2 minutes to calculate
    %     BestGroupOFA=1:4;
    BestGroupOFA=BestGroup;
    [RelaxOFA{1:3}]=CalcT1byFAfw2(FFA2DData(:,BestGroupOFA)',FAs(BestGroupOFA),TRs(BestGroupOFA));
    FRelax4D3=Reshape2DCto4D(RelaxOFA,FBrainMaskS>-1);
    
    Raw2Nii(squeeze(FRelax4D3(:,:,:,1)),[RWorkingP 'T13DOFA.nii'],'float32', MeanFN);
    Raw2Nii(squeeze(FRelax4D3(:,:,:,2)),[RWorkingP 'PD3DOFA.nii'],'float32', MeanFN);
    Raw2Nii(squeeze(FRelax4D3(:,:,:,3)),[RWorkingP 'RMS3DOFA.nii'],'float32', MeanFN);
    
    MidSli=floor(size(FA4D,3)/2);
    load(PrepareFN,'BadSlicesF2');
    GoodSlices=setdiff(1:size(FA4D,3),BadSlicesF2);
    
    % ASK GILAD - What is the purpose of the following paragraph?
    %             Currently the GUI does not define that parameter.
    % ASNWER -  Correction for flip angles (exp. � machine says its 20 degrees while the real one is 20.2).
    %                     Implemented under CalcRelaxForVolNFA.m. Currently not used.
    ExtractFAs=Options.ExtractFAs && numel(BestGroup)>2;
    if(ExtractFAs)
        BadNFAs=@(x) any(x==0.5) || sum((x-floor(x))==0)>1;
        if(exist(MatFN,'file') && ~CalcForce)
            disp('Skipping FA estimation');
        else
            disp(WorkingP);
            Tmp=FA4D;
            for ii=1:nFAs
                Tmp(:,:,:,ii)=FA4D(:,:,:,ii)/(T2SFactors(ii).*GainFactors(ii));
            end
            ONFAsx=CalcRelaxForVolNFAf2(Tmp(:,:,:,BestGroup),FBrainMaskS,FAs(BestGroup),TRs(BestGroup),10,1000);
            ONFAs=FAs;
            ONFAs(BestGroup)=ONFAsx;
            ONFAsBefore=ONFAs;
            if(BadNFAs(ONFAsx))
                disp('Extracted bad FAs!');
                AddToLog(WorkingP,'b_5badfa',['Extracted bad FAs! ' num2str(ONFAs)],[],4);
                %                 Error('Extracted bad FAs!');
                ONFAs=FAs;
            end
        end
    end
    
    [RelaxNFA{1:3}]=CalcT1byFAfw2(FFA2DData(:,BestGroup)',ONFAs(BestGroup),TRs(BestGroup));
    FRelax4DN=Reshape2DCto4D(RelaxNFA,FBrainMaskS>-1);
    
    if(ExtractFAs)
        NotInBestGroup=setdiff(1:nFAs,BestGroup);
        for ii=1:numel(NotInBestGroup)
            CurFA=NotInBestGroup(ii);
            ONFAs(CurFA)=CalcFAgSigT1M0(FFA2DData(FFBrainMaskSF,CurFA),RelaxNFA{1}(FFBrainMaskSF),RelaxNFA{2}(FFBrainMaskSF),TRsF(CurFA),100,max(0.5,FAsF(CurFA)/2),FAsF(CurFA)*2);
        end
    end
    
    AddToLog(WorkingP,'b_3p1',['ONFAsBefore:  ' num2str(ONFAsBefore,'% 8.2f')]);
    AddToLog(WorkingP,'b_3p2',['ONFAs:        ' num2str(ONFAs,'% 8.2f')]);
    %%
    %         FMsk=find(FBrainMaskS);
    %         CostFunc=@(x) gCost(FFA2DData(FMsk,end),x(2)*SPGRfM(RelaxNFA{1}(FMsk)',RelaxNFA{2}(FMsk)',x(1),TRs(end)),'RMS');
%             CostFunc=@(x) gCost(FFA2DDataX(FMsk,end),SPGRfM(RelaxNFA{1}(FMsk)',RelaxNFA{2}(FMsk)',x(1),TRs(end)),'RMS');
%             BestFA2=patternsearch(CostFunc,[FAs(end)])
%             ONFAs(end)=BestFA2

    %         BestFAFac=BFGSf(CostFunc,[FAs(end) 1],[0.5 0.5],[FAs(end)*3 2],struct('GradObj','off','Display','off'))
    %         BestFAFac2=patternsearch(CostFunc,[FAs(end) 1])
    %
    %         ONFAs=[ONFAsBefore BestFAFac(1)];
    %         Factorsx=Factors;
    %         Factorsx(end)=BestFAFac(2);
    %% EM try - seems good
    EMtry=false;
    if(EMtry)
        FFA2DDataX=Reshape4d22d(FA4D,FBrainMaskS>-1);
        FFA2DDataX=repMulti(FFA2DDataX,1./T2SFactors);
        
        BestGroupNewFAs=BestGroup;
        for ii=1:3
            RelaxNewFA{ii}=RelaxOFA{ii}(FFBrainMaskSF);
        end
        disp(median(RelaxNewFA{3}(:)));
        disp(median(RelaxOFA{3}(FFBrainMaskSF)));
        disp(median(RelaxNFA{3}(FFBrainMaskSF)));
        for tt=1:5
            for ii=1:nFAs
                NewFAs(ii)=CalcFAgSigT1M0(FFA2DDataX(FFBrainMaskSF,ii),RelaxNewFA{1}(:),RelaxNewFA{2}(:),TRsF(ii),100,max(0.5,FAsF(ii)/2),FAsF(ii)*2);
            end
            [RelaxNewFA{1:3}]=CalcT1byFAfw2(FFA2DDataX(FFBrainMaskSF,BestGroupNewFAs)',NewFAs(BestGroupNewFAs),TRs(BestGroupNewFAs));
            disp(median(RelaxNewFA{3}(:)))
        end
    end
    %%
    figure(8181);clf;
    for i=1:nFAs
        gsubplot(5,nFAs,1,i);
        Tmp=squeeze(FA4D(:,:,:,i));
        [Mn Mx]=FindDR(Tmp(FBrainMaskS));
        imagesc(mritransform(Tmp(:,:,MidSli)),[0 Mx]);
        title(['[0 ' num2str(Mx) ']']);
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        ylabel([num2str([FAs(i) ONFAs(i)])]);
        xlabel(num2str([TRs(i) TEs(i)]));
        title('RAW Data');
        gsubplot(5,nFAs,2,i);
        Tmp(:)=SPGRfM(RelaxOFA{1}',RelaxOFA{2}',FAs(i),TRs(i))*T2SFactors(i)*GainFactors(i);
%         Tmp(:)=SPGRfM(RelaxOFA{1}',RelaxOFA{2}',FAs(i),TRs(i))*T2SFactors(i)*Gains3(i);
%         Tmp(:)=SPGRfM(RelaxOFA{1}',RelaxOFA{2}',FAs(i),TRs(i))*T2SFactors(i);
        imagesc(mritransform(Tmp(:,:,MidSli)),[0 Mx]);
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        title('Sim OFA');
        gsubplot(5,nFAs,3,i);
        Tmp=squeeze(FA4D(:,:,:,i))-Tmp;
        imagesc(mritransform(Tmp(:,:,MidSli)),[-100 100]);
        title(['Diff OFA [-100 100]']);
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        xlabel(median(abs(Tmp(FBrainMaskS))));
        gsubplot(5,nFAs,4,i);
        %             Tmp(:)=SPGRfM(RelaxNFA{1}',RelaxNFA{2}',ONFAs(i),TRs(i))*T2SFactors(i)*GainFactors(i);
        Tmp(:)=SPGRfM(RelaxNFA{1}',RelaxNFA{2}',ONFAs(i),TRs(i))*T2SFactors(i);
        %             Tmp(FBrainMaskS)=SPGRfM(RelaxNewFA{1}',RelaxNewFA{2}',NewFAs(i),TRs(i))*T2SFactors(i);
        imagesc(mritransform(Tmp(:,:,MidSli)),[0 Mx]);
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        title('Sim NFA');
        gsubplot(5,nFAs,5,i);
        Tmp=squeeze(FA4D(:,:,:,i))-Tmp;
        imagesc(mritransform(Tmp(:,:,MidSli)),[-100 100]);
        title(['Diff NFA [-100 100]']);
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        xlabel(median(abs(Tmp(FBrainMaskS))));
    end
    %%
    saveas(8181,[WorkingP 'T1andFAestimationSlice.png']);
    saveas(8181,[WorkingP 'T1andFAestimationSlice.fig']);
    close(8181);
    AddToLog(WorkingP,'yb_5b',['Img T1 fit for Mid slice, Extracted FAs: ' num2str(ONFAs)],'T1andFAestimationSlice.png');
    %%
    Raw2Nii(squeeze(FRelax4DN(:,:,:,1)),[RWorkingP 'T13DNFA.nii'],'float32', MeanFN);
    Raw2Nii(squeeze(FRelax4DN(:,:,:,2)),[RWorkingP 'PD3DNFA.nii'],'float32', MeanFN);
    Raw2Nii(squeeze(FRelax4DN(:,:,:,3)),[RWorkingP 'RMS3DNFA.nii'],'float32', MeanFN);
    
    T1Res{CurCGroup,2,1}=[RWorkingP 'T13DNFA.nii'];
    
    OFASli=squeeze(FRelax4D3(:,:,MidSli,1));
    NFASli=squeeze(FRelax4DN(:,:,MidSli,1));
    OFARMS=squeeze(FRelax4D3(:,:,MidSli,3));
    NFARMS=squeeze(FRelax4DN(:,:,MidSli,3));
    MaxT1x=max([max(OFASli(:)) max(NFASli(:))]);
    MaxT1x=3000;
    MaxRMSx=max([max(OFARMS(:)) max(NFARMS(:))]);
    [MnV MaxRMSx]=FindDR(OFARMS(OFARMS>0));
    figure(91919191);clf;
    subplot(2,2,1);imagesc(mritransform(OFASli),[0 MaxT1x]);title('T_1 nominal FAs'); xlabel(num2str(FAsF));
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
    subplot(2,2,2);imagesc(mritransform(NFASli),[0 MaxT1x]);title('T_1 changed FAs'); xlabel(num2str(ONFAs,'% 2.2f'));
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
    subplot(2,2,3);imagesc(mritransform(log(OFARMS)),[0 log(MaxRMSx)]);title('log RMS nominal FAs');
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
    subplot(2,2,4);imagesc(mritransform(log(NFARMS)),[0 log(MaxRMSx)]);title('log RMS extracted FAs');
    set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
    xlabel(num2str(log(MaxRMSx)));
    saveas(91919191,[WorkingP 'T1andFAestimation.png']);
    saveas(91919191,[WorkingP 'T1andFAestimation.fig']);
    close(91919191);
    save(MatFN,'ONFAs','TRsF','FAsF','ONFAsBefore');
    EEBrainMask=imerode(EBrainMask,se);
    EEEBrainMask=imerode(EEBrainMask,se);
    AddToLog(WorkingP,'yb_5c',['Img T1 fit RMS'],'T1andFAestimation.png');
    %%
    ForSeg=loadniidata(T1Res{CurCGroup,2,1});
    ForSeg(~BrainMask)=NaN;
    ForSegFN=[RWorkingP 'ForSeg.nii'];
    Raw2Nii(ForSeg,ForSegFN,'float32', MeanFN);
    
    T1NFA_SegP=SPM_Segment(ForSegFN,false,[],false,FMaskFN);
    %
    %         FMaskFN2=[RWorkingP 'MaskForSeg.nii'];
    %         BrainMaskX=BrainMask*0+1;
    %         BrainMaskX(~BrainMask)=NaN;
    %         Raw2Nii(BrainMask,FMaskFN2,'uint8', MeanFN);
    %
    %         T1NFA_SegP=SPM_Segment(T1Res{CurCGroup,2,1},true,[],false,FMaskFN2);
    
    %         T1NFA_SegP=SPM_Segment(T1Res{CurCGroup,2,1},false,[],false,FMaskFN);
    
    %         T1NFA_SegP=SPM_Segment(T1Res{CurCGroup,2,1},true,[],false);
    
    T1Seg3D(:,:,:,1)=loadniidata([T1NFA_SegP 'c1ForSeg.nii'])/256;
    T1Seg3D(:,:,:,2)=loadniidata([T1NFA_SegP 'c2ForSeg.nii'])/256;
    T1Seg3D(:,:,:,3)=loadniidata([T1NFA_SegP 'c3ForSeg.nii'])/256;
    T1Seg3D(:,:,:,4)=loadniidata(FMaskFN);
    [Tmp, T1Seg3DAll]=max(T1Seg3D(:,:,:,1:3),[],4);
    T1Seg3DAll(~T1Seg3D(:,:,:,4))=0;
    T1Cleaned=loadniidata([T1NFA_SegP 'mForSeg.nii']);
    
    CSFMask=T1Seg3DAll==3 & T1Seg3D(:,:,:,3)>0.99 & EEEBrainMask==1;
    WMMask=T1Seg3DAll==2 & T1Seg3D(:,:,:,2)>0.99 & EEEBrainMask==1;
    
    T1Seg3DAllx=T1Seg3DAll;
    T1Seg3DAllx(CSFMask)=4;
    T1Seg3DAllx(WMMask)=5;
    
    Tmp=max(FBrainMaskS,[],3);
    F=find(max(Tmp,[],2));
    GoodRows=F(1):F(end);
    F=find(max(Tmp,[],1));
    GoodCols=F(1):F(end);
    for i=1:numel(GoodSlices)
        CurSli=GoodSlices(i);
        I=squeeze(T1Cleaned(:,:,CurSli));
        %         [Tmp MaxV]=FindDR(I);
        MaxV=4000;
        ClrM=[1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1];
        IRGB=repmat(min(1,I/MaxV),[1 1 3]);
        for tt=1:3
            BW2 = bwmorph(squeeze(T1Seg3DAll(:,:,CurSli)==tt),'remove');
            for kk=1:3
                TmpI=IRGB(:,:,kk);
                TmpI(BW2)=ClrM(tt,kk);
                IRGB(:,:,kk)=TmpI;
            end
        end
        for tt=4:5
            BW2 = bwmorph(squeeze(T1Seg3DAllx(:,:,CurSli)==tt),'remove');
            for kk=1:3
                TmpI=IRGB(:,:,kk);
                TmpI(BW2)=ClrM(tt,kk);
                IRGB(:,:,kk)=TmpI;
            end
        end
        IRGB3(:,:,:,i)=IRGB;
    end
    figure(9899);clf;
    montage(mritransform(IRGB3(GoodRows,GoodCols,:,:)))
    title(num2str(GoodSlices));
    xlabel(num2str([FAsF; ONFAs],'% 2.2f'));
    saveas(9899,[WorkingP 'T1Seg.png']);
    saveas(9899,[WorkingP 'T1Seg.fig']);
    close(9899);
    AddToLog(WorkingP,'yb_6',['Img segmentation. Red - GM, Green - WM, Blue - CSF, Magenta - WM for reference, Yellow - CSF for reference.'],'T1Seg.png');
    %%
    Raw2Nii(CSFMask,[WorkingP 'RefAuto_CSF_2430.nii'],'float32', MeanFN);
    Raw2Nii(WMMask,[WorkingP 'RefAuto_WM_830.nii'],'float32', MeanFN);
end

a=whos;
BigMem={a([a.bytes]>10000).name}';
clear(BigMem{:});
AddToLog(WorkingP,'b_9','Relaxometry finished');
save(RelaxFN);
disp(['DCET1_RelaxForSubjectf finished for ' WorkingP]);