function DCET1_RelaxForSubjectf(DCEInfos,WorkingP,DoN3,DoGlobal,DoDeviations,CalcForce,WhichMean,Options,MainInfo)

LogFN=[WorkingP 'Log.mat'];
if(~exist(LogFN,'file'))
    SN=[WorkingP];
    Log.a_00={['\\title{' SN '}\r\n\\maketitle\r\n']};
    save(LogFN,'Log');
end

JustT1Mode=false;
if(isempty(MainInfo))
    JustT1Mode=true;
end

AddToLog(WorkingP,'b_00','\\subsection*{Relaxometry}');

% DoN3,DoGlobal,DoDeviations
RelaxFN=[WorkingP 'Relax.mat'];
%% Create the Relaxometry directory and needed initial structs in case it does not exist already
DCERelaxP=[WorkingP 'Relaxometry' filesep];
if(~isempty(MainInfo))
    DCEInfos=DCEInfos(~ismember({DCEInfos.SeriesInstanceUID},MainInfo.SeriesInstanceUID));
end

% Create relaxometry directory
mkdir(DCERelaxP);
%% Get the number of flip angle sets
% Create a copy of DCE Infos from which we pick every iteration the relevant infos
Copy_DCEInfos = DCEInfos;
Tmp=[Copy_DCEInfos.ImagesInAcquisition];
Philips=(Tmp(end)==25);
Tmp(Tmp==25)=2500;
SingleSPGRIdxs=find(Tmp>1.5*(Copy_DCEInfos(1).ImagesInAcquisition));
%   NOTE - I assume all the sets have the same angle scans (repeated for each set)

% First angle used
first_angle = DCEInfos(1).FlipAngle;
% All indices with same angle
index_new_sets = find([DCEInfos.FlipAngle]==first_angle);
index_new_sets = union(index_new_sets,SingleSPGRIdxs);
% Number of sets
num_FAs_sets = numel(index_new_sets)-Philips;
StartIdxs=index_new_sets;
EndIdxs=[index_new_sets(2:end)-1 size(Copy_DCEInfos,2)];

FAsFF=[DCEInfos.FlipAngle];
TEsFF=[DCEInfos.EchoTime];
TRsFF=[DCEInfos.RepetitionTime];
R1sFF=[DCEInfos.R1];
R2sFF=[DCEInfos.R2];
TGsFF=[DCEInfos.TG];

FAsBySet=FAsFF(StartIdxs);
TEsBySet=TEsFF(StartIdxs);
TRsBySet=TRsFF(StartIdxs);
R1sBySet=R1sFF(StartIdxs);
R2sBySet=R2sFF(StartIdxs);
TGsBySet=[];
if(numel(TGsFF)>=max(StartIdxs))
    TGsBySet=TGsFF(StartIdxs);
end

AddToLog(WorkingP,'b_00001',['FAsBySet: ' num2str(FAsBySet)]);
AddToLog(WorkingP,'b_00002',['TEsBySet: ' num2str(TEsBySet)]);
AddToLog(WorkingP,'b_00003',['TRsBySet: ' num2str(TRsBySet)]);
AddToLog(WorkingP,'b_00004',['R1sBySet: ' num2str(R1sBySet)]);
AddToLog(WorkingP,'b_00005',['R2sBySet: ' num2str(R2sBySet)]);
AddToLog(WorkingP,'b_00006',['TGsBySet: ' num2str(TGsBySet)]);

%% Repeat the T1 process for every set of angles
FASetsToWorkOn=1:num_FAs_sets;

for j =FASetsToWorkOn
    
    % Create a directory for each set
    Set_Dir =[DCERelaxP 'Set_Num_' sprintf('%d',j) filesep];
    mkdir(Set_Dir);
    % Set the working directory to be the same
    RWorkingP=Set_Dir;
    
    % Get the relevant DCE infos to work on
    first_angle_to_work = StartIdxs(j);
    last_angle_to_work   = EndIdxs(j);
    
    % If we have less than 2 angles to work with, we can not create T1
    % map -> Error!
    
    
    % Handle single Angle!
    if (last_angle_to_work - first_angle_to_work + 1 < 2)
        
        DCEInfos = Copy_DCEInfos(first_angle_to_work :last_angle_to_work);
        % Get average Time Stamp of every scan
        Time_Stamps =   {DCEInfos.SeriesDateTime};
        size_of_list = size(Time_Stamps,2);
        % Get the middle index to ve the time stamp
        middle_index = ceil(size_of_list/2);
        Mid_Time_Stamp = [Time_Stamps(middle_index)];
        
        % Create a text file which will hold the average time of the scan
        file_name = [RWorkingP 'Middle_Time_Stamp.mat'];
        save (file_name,'Mid_Time_Stamp')
        
        % Creating nifti for all original FAs
        % Create file name for nifti
        CurOtherFAsInfos=DCEInfos;
        i=1;
        CurFA=round(CurOtherFAsInfos(i).FlipAngle);
        OtherFABaseFNs{i}=[Set_Dir 'OrigFA_' sprintf('%02d',CurFA) '_' sprintf('%02d',CurOtherFAsInfos(i).SeriesNumber) '.nii'];
        
        % Create temp nifti out of dicom data
        Tmp=gDicom2Nifti(CurOtherFAsInfos(i).Path,OtherFABaseFNs{i});
        
        % Same nifti file name concatenated with "_0001"
        MainSubFN=[Set_Dir 'OrigFA_' sprintf('%02d',CurFA) '_' sprintf('%02d',CurOtherFAsInfos(i).SeriesNumber) '_0001.nii'];
        
        % Not clear what the purpose of the following...
        if(Tmp>1 && Tmp<10 && exist(MainSubFN,'file'))
            MainSubFN=[Set_Dir 'OrigFA_' sprintf('%02d',CurFA) '_' sprintf('%02d',CurOtherFAsInfos(i).SeriesNumber) '_' sprintf('%04d',floor(Tmp/2))   '.nii'];
            movefile(MainSubFN,OtherFABaseFNs{i});
            delete([Set_Dir 'OrigFA_' sprintf('%02d',CurFA) '_' sprintf('%02d',CurOtherFAsInfos(i).SeriesNumber) '_*.nii']);
        end
        
        if(~exist(OtherFABaseFNs{i},'file'))
            error('DCET1_RelaxForSubjectf Dcm2Nii NotCreated');
        end
        
        %% Coreg dest
        CMeanFN=[WorkingP 'CopyDCEMean.nii'];
        MeanFN=[WorkingP 'DCEMean.nii'];
        delete(CMeanFN);
        try
            copyfile(MeanFN,CMeanFN,'f');
        catch MErr
        end
        
        if(~exist(CMeanFN,'file'))
            error('DCET1_RelaxForSubjectf: Problem copying CMeanFN');
        end
        
        %% Coreg
        
        OtherFAMatFNs{i}=CoregEstimate(OtherFABaseFNs{i},CMeanFN,CalcForce);
        OtherFACrgFNs{i}=CoregWrite(OtherFABaseFNs{i},OtherFAMatFNs{i},CalcForce,Set_Dir,false,CMeanFN);
        delete(CMeanFN);
        
        %% ---- End single angle ------
        
        
        %         helpdlg('Error in Data!!! One of the DESPOT T1 data has less than 2 flip angles. Can not create T1 map this way!');
        AddToLog(WorkingP,['b_00a' num2str(j)],['Skipping DESPOT1 set ' num2str(j) ' vols ' num2str(last_angle_to_work) ':' num2str(first_angle_to_work)]);
        %             error('Error in Data!!! One of the DESPOT T1 data has less than 2 flip angles. Can not create T1 map this way!');
        continue;
    end
    
    DCEInfos = Copy_DCEInfos(first_angle_to_work :last_angle_to_work);
    
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
        CurFA=round(CurOtherFAsInfos(i).FlipAngle);
        
        % Create file name for nifti
        OtherFABaseFNs{i}=[Set_Dir 'OrigFA_' sprintf('%02d',CurFA) '_' sprintf('%02d',CurOtherFAsInfos(i).SeriesNumber) '.nii'];
        
        % Create temp nifti out of dicom data
        Tmp=gDicom2Nifti(CurOtherFAsInfos(i).Path,OtherFABaseFNs{i});
        % Same nifti file name concatenated with "_01"
        MainSubFN=[Set_Dir 'OrigFA_' sprintf('%02d',CurFA) '_' sprintf('%02d',CurOtherFAsInfos(i).SeriesNumber) '_01.nii'];
        if(~exist(MainSubFN,'file'))
            MainSubFN=[Set_Dir 'OrigFA_' sprintf('%02d',CurFA) '_' sprintf('%02d',CurOtherFAsInfos(i).SeriesNumber) '_0001.nii'];
        end
        if(Tmp==3)
            MainSubFN=[Set_Dir 'OrigFA_' sprintf('%02d',CurFA) '_' sprintf('%02d',CurOtherFAsInfos(i).SeriesNumber) '_0002.nii'];
        end
        
        % Not clear what the purpose of the following...
        if(Tmp>1 && Tmp<10 && exist(MainSubFN,'file'))
            
            movefile(MainSubFN,OtherFABaseFNs{i});
            delete([Set_Dir 'OrigFA_' sprintf('%02d',CurFA) '_' sprintf('%02d',CurOtherFAsInfos(i).SeriesNumber) '_*.nii']);
            
        end
        
        if(~exist(OtherFABaseFNs{i},'file'))
            error('DCET1_RelaxForSubjectf Dcm2Nii NotCreated');
        end
        
    end
    
    %% Coreg dest
    
    MedFA=median(FAsF);
    [tmp, MedFAI]=min(abs(FAsF-MedFA));
    
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
    
    if(~strcmp(WhichMean,'No coreg'))
        if(~exist(CMeanFN,'file'))
            error('DCET1_RelaxForSubjectf: Problem copying CMeanFN');
        end
    end
    
    %% Coreg
    if(strcmp(WhichMean,'No coreg'))
        MeanFN=[WorkingP 'DCEMean.nii'];
        for i=1:nOtherFAs
            OtherFACrgFNs{i}=OtherFABaseFNs{i};
        end
    else
        for i=1:nOtherFAs
            OtherFAMatFNs{i}=CoregEstimate(OtherFABaseFNs{i},CMeanFN,CalcForce);
            OtherFACrgFNs{i}=CoregWrite(OtherFABaseFNs{i},OtherFAMatFNs{i},CalcForce,Set_Dir,false,CMeanFN);
        end
    end
    delete(CMeanFN);
    disp('Relaxometry coregistration finished');
    
    %% Full brain mask for relaxometry
    PrepareFN=[WorkingP 'AfterPrepare4D.mat'];
    
    %load(PrepareFN,'BrainMask','C2','RepSli','MeanVol');
    if(JustT1Mode)
        MeanVol=loadniidata(MeanFN);
        Raw2Nii(MeanVol,CMeanFN,'float32',MeanFN);
        if(Options.Mask_thresh<0)
            [BetStrippedMeanFN, BetMaskFN]=mcbet2(CMeanFN,CalcForce,abs(Options.Mask_thresh));
        else
            [BetStrippedMeanFN, BetMaskFN]=mcbet2(CMeanFN,CalcForce);
        end
        ManMaskFN=[WorkingP 'Manual_BrainMask.nii'];
        if(exist(ManMaskFN,'file'))
            %AddToLog(WorkingP,'a_2caaaaa','Using Manual_BrainMask');
            AddToLog(WorkingP,'a_2caaaaa','Using Manual.BrainMask');
            BrainMask=loadniidata(ManMaskFN)>0;
        else
            BrainMask=loadniidata(BetMaskFN)>0;
        end
    else
        load(PrepareFN,'BrainMask');
    end
    
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
    AddToLog(WorkingP,['b_2d' num2str(j)],[num2str(j) ' TRs:        ' num2str(TRsF,'% 8.2f')]);
    AddToLog(WorkingP,['b_2e' num2str(j)],[num2str(j) ' TEs:        ' num2str(TEsF,'% 8.2f')]);
    AddToLog(WorkingP,['b_2f' num2str(j)],[num2str(j) ' FAs:        ' num2str(FAsF,'% 8.2d')]);
    AddToLog(WorkingP,['b_2g' num2str(j)],[num2str(j) ' R1s:        ' num2str(R1sF,'% 8.2d')]);
    AddToLog(WorkingP,['b_2h' num2str(j)],[num2str(j) ' R2s:        ' num2str(R2sF,'% 8.2d')]);
    AddToLog(WorkingP,['b_2i' num2str(j)],[num2str(j) ' TGs:        ' num2str(TGsF,'% 8.2d')]);
    
    T2SApproximation=16;
    AddToLog(WorkingP,['b_2j' num2str(j)],[num2str(j) ' Using $T_2^*$ approximation of ' num2str(T2SApproximation) 'ms']);
    T2SFactors=exp(-TEsF/T2SApproximation);
    AddToLog(WorkingP,['b_2k' num2str(j)],[num2str(j) ' $T_2^*$ Factors: ' num2str(T2SFactors,'% 8.2f')]);
    T2SSigPerDif=100*(sqrt(max(T2SFactors)/min(T2SFactors))-1);
    AddToLog(WorkingP,['b_2l' num2str(j)],[' Estimated signal change because of echo time differences ' num2str(T2SSigPerDif) '\\%']);
    %%
    TRT1difF=100*(sqrt(max(TRsF)/min(TRsF))-1);
    AddToLog(WorkingP,['b_2m' num2str(j)],[num2str(j) ' Estimated $T_1$ inaccuracy because of repetition time differences on all ' num2str(TRT1difF,'%2.2f') '\\%']);
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
    AddToLog(WorkingP,['b_2n' num2str(j)],[num2str(j) ' Chosen group to for FA extraction: ' num2str(BestGroup,'% d')]);
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
    warning('off','MATLAB:illConditionedMatrix');
    warning('off','MATLAB:singularMatrix');
    warning('off','MATLAB:nearlySingularMatrix');
    
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
    
    TRGroup= TRGroupC{CurCGroup};
    Needed = TRGroup; % [1:4];
    FAs    = FAsF(Needed);
    TRs    = TRsF(Needed);
    TEs    = TEsF(Needed);
    FA4D   = FA4DF(:,:,:,Needed);
    nFAs   = numel(FAs);
    SFA    = size(FA4D);
    ANFAs  = zeros(SFA(3),nFAs);
    ValMsk = all(FA4D>100,4);
    % FBrainMaskS=FBrainMask & ValMsk;
    
    % CalcRelaxForVolNFA;
    if(~exist('ONFAs','var'))
        ONFAs=FAs;
    end
    %     load(MatFN,'ONFAs','FFA2DData','FBrainMaskS');
    %     load(MatFN,'ONFAs');
    
    FBrainMaskS=BrainMask;
    
    FFA2DData = Reshape4d22d(FA4D,FBrainMaskS>-1);
    FFA2DData = repMulti(FFA2DData,1./T2SFactors);
    
    %% Gain factors -
    if (Options.Calc_Gains_Diff)
        %  Calculate the signal gain for each volume of the DESPOT base FAs
        
        % Initiate all gains with '1'
        GainFactors   = ones(1,nFAs);
        % Take linear indices of the brain mask
        FFBrainMaskSF = find(FBrainMaskS);
        % Take 1000 random pixels from the brain volume
        Tmp           = FFA2DData(FFBrainMaskSF,:);
        Tmp           = Tmp(getKrandomSamples(size(Tmp,1),1000),:);
        clear AllEstGainFactors EstFac AllGainFactors
        % Initiate iterations parameters
        num_iter          = 4;
        EstFac            = zeros(1,nFAs);
        AllEstGainFactors = zeros(num_iter,nFAs);
        AllGainFactors    = zeros(num_iter,nFAs);
        
        % After 2 iterations it should converge
        for iter_idx = 1:num_iter
            
            % For each angle, calculate a simulated signal out of the rest of
            % the angles
            for ii=1:nFAs
                % Take the rest of the angles
                Others         = setdiff(1:nFAs,ii);
                
                % Get the set number we are working on
                Set_Num = j;
                
                % Regular DESPOT calculation when this is the first set or
                % when we don't want to use a single M0
                if (Set_Num == 1 || ~Options.Use_Single_M0)
                    
                    % Calculate T1, PD, RMS for by rest of angles
                    [RelaxCG{1:3}] = CalcT1byFAfw2(Tmp(:,Others)',FAs(Others),TRs(Others));
                    
                    % Save the first M0 calculated
                    First_M0_Map_CG = RelaxCG{2};
                else
                    
                    % User wanted single M0 for all additional DESPOT calculations
                    % Make sure there is only 1 angle in that despot
                    if ( numel(FAs) == 1 )
                        [RelaxCG{1:3}] = CalcT1byFAfw2_Single_Angle(Tmp(:)',FAs,TRs,First_M0_Map_CG);
                    else
                        error('-E- There is more than 1 angle in DESPOT when trying to calculate T1');
                    end
                    
                end
                
                % Create the simulated signal
                SimByNomFA     = SPGRfM(RelaxCG{1}',RelaxCG{2}',FAs(ii),TRs(ii));
                % Take the median gain value between the original signal and
                % simulated one
                EstFac(ii)     = median(Tmp(:,ii)./SimByNomFA);
            end
            
            % Put the gain in the relevent iteration index
            AllEstGainFactors(iter_idx,:) = EstFac;
            % Round the 2^x numbers
            CurGainFactors                = 2.^round(log2(EstFac));
            % Save the rounded gain value in the relevant iteration index
            AllGainFactors(iter_idx,:)    = CurGainFactors;
            % Calculate a new value for the signal (normalized by the gain
            Tmp                           = repMulti(Tmp,1./CurGainFactors);
        end
        NoUseCurrently_GainFactors=prod(AllGainFactors,1);
        Gains3=NoUseCurrently_GainFactors.*AllEstGainFactors(end,:);
        AddToLog(WorkingP,['b_3d' num2str(j)],[num2str(j) ' Estimated gains, not used but check:        ' num2str(AllEstGainFactors(1,:),'% 8.2f')]);
        AddToLog(WorkingP,['b_3e' num2str(j)],[num2str(j) ' Rounded gains, not used but check:        ' num2str(NoUseCurrently_GainFactors,'% 8.2f')],[],2);
        
        % Each angle gain is the product of all iterations rounded gains
        NoUseCurrently_GainFactors = prod(AllGainFactors,1);
        %Gains3                     = NoUseCurrently_GainFactors.*AllEstGainFactors(end,:);
        
        % Add to log file the estimated and rounded gains
        AddToLog(WorkingP,'b_3n1',['Estimated gains, not used but check:      ' num2str(AllEstGainFactors(1,:),'% 8.2f')]);
        AddToLog(WorkingP,'b_3n2',['Rounded gains, not used but check:        ' num2str(NoUseCurrently_GainFactors,'% 8.2f')],[],2);
        
    end
    
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
    
    % Get the set number we are working on
    Set_Num = j;
    
    % Regular DESPOT calculation when this is the first set or
    % when we don't want to use a single M0
    if (Set_Num == 1 || ~Options.Use_Single_M0)
        
        % Calculate T1, PD, RMS for by rest of angles
        [RelaxOFA{1:3}] = CalcT1byFAfw2(FFA2DData(:,BestGroupOFA)',FAs(BestGroupOFA),TRs(BestGroupOFA));
        
        % Save the first M0 calculated
        First_M0_Map_OFA = RelaxOFA{2};
        
    else
        % User wanted single M0 for all additional DESPOT calculations
        % Make sure there is only 1 angle in that despot
        if ( numel(FAs) == 1 )
            [RelaxOFA{1:3}]=CalcT1byFAfw2_Single_Angle(FFA2DData(:,BestGroupOFA)',FAs(BestGroupOFA),TRs(BestGroupOFA),First_M0_Map_OFA);
        else
            error('-E- There is more than 1 angle in DESPOT when trying to calculate T1');
        end
    end
    
    FRelax4D3=Reshape2DCto4D(RelaxOFA,FBrainMaskS>-1);
    %     MidSli=floor(size(FA4D,3)/2);
    %     figure;imagesc(mritransform(FRelax4D3(:,:,MidSli,1)),[0 4000]);
    
    Raw2Nii(squeeze(FRelax4D3(:,:,:,1)),[RWorkingP 'T13DOFA.nii'],'float32', MeanFN);
    Raw2Nii(squeeze(FRelax4D3(:,:,:,2)),[RWorkingP 'PD3DOFA.nii'],'float32', MeanFN);
    Raw2Nii(squeeze(FRelax4D3(:,:,:,3)),[RWorkingP 'RMS3DOFA.nii'],'float32', MeanFN);
    
    MidSli=floor(size(FA4D,3)/2);
    if(JustT1Mode)
        SDCE=size(FA4D);
        %  Mark the first and the last slices to be false
        BadSlicesA=[1 size(MeanVol,3)];
        
        % Remove the first and last slices
        BrainMask(:,:,BadSlicesA)=false;
        
        %Selecting slices
        disp('Find problematic slices');
        
        MedSlice=zeros(SDCE(3),1);
        % Go over every slice
        for i=1:SDCE(3)
            % Slice mask at every iteration holds values different than 0 for the
            % current slice only (this is why we initiate it with zeros)
            SliceMsk=zeros(SDCE(1),SDCE(2),SDCE(3));
            SliceMsk(:,:,i)=BrainMask(:,:,i);
            
            % Take the first time index only because we want an image before the enhancement.
            % We should have taken the entire baseline (all images before
            % enhancmenet), but this is an approximation
            Tmp=Reshape4d22d(MeanVol,SliceMsk);
            % Get the median value of every slice
            MedSlice(i,:)=median(Tmp,1);
        end
        
        % The following function is problematic (dimenstions for repPlus,repMulti
        % are incorrect)
        %NMedSlice=repMulti(repPlus(MedSlice,-mean(MedSlice,1)),1./std(MedSlice,0,1));
        
        % BadSlices=any(abs(NMedSlice)>2,2);
        NaNSlices=isnan(MedSlice);
        disp('Looking for weird slices');
        if(numel(MedSlice(~NaNSlices))<2)
            error('Only one slice!');
        end
        
        BadSlices=NaNSlices;
        
        % The following function is checking if there are any other bad slices
        % A slice is considered bad if its median is far (>50) than the 3 middle slices (which should have a reliable signal
        BadSlices(~NaNSlices)=abs(MedSlice(~NaNSlices)-mean(MedSlice(MidSli-1:MidSli+1)))>50-Philips*30;
        BadSlicesF=find(BadSlices);
        disp(['Ignoring slices ' num2str(BadSlicesF')]);
        AddToLog(WorkingP,'a_2ccc',['Ignoring slices ' num2str(BadSlicesF')]);
        figure(78362);clf;subplot(1,2,1);
        % Plot the median value of each slice
        plot(1:SDCE(3),MedSlice,'b',BadSlicesF,MedSlice(BadSlicesF),'ro');
        title('Median of slices and ignored ones');
        xlabel('Slices');
        ylabel('Median');
        MskSlices=ones(size(MeanVol))>0;
        MskSlices(:,:,BadSlicesF)=false;
        % Get a list of all bad slices (unite the extreme slices and the ones >50)
        BadSlicesF2=union(BadSlicesF,[1 size(MeanVol,3)]);
    else
        load(PrepareFN,'BadSlicesF2');
    end
    
    GoodSlices=setdiff(1:size(FA4D,3),BadSlicesF2);
    
    % ASK GILAD - What is the purpose of the following paragraph?
    %             Currently the GUI does not define that parameter.
    % ASNWER -  Correction for flip angles (exp. MRI machine says its 20 degrees while the real one is 20.2).
    %                     Implemented under CalcRelaxForVolNFA.m. Currently not used.
    ExtractFAs=Options.ExtractFAs && numel(BestGroup)>2 && ~Philips;
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
                AddToLog(WorkingP,['b_5badfa' num2str(j)],[num2str(j) ' Extracted bad FAs! ' num2str(ONFAs)],[],4);
                %                 Error('Extracted bad FAs!');
                ONFAs=FAs;
            end
        end
    end
    
    % Get the set number we are working on
    Set_Num = j;
    
    % Regular DESPOT calculation when this is the first set or
    % when we don't want to use a single M0
    if (Set_Num == 1 || ~Options.Use_Single_M0)
        
        % Calculate T1, PD, RMS for by rest of angles
        [RelaxNFA{1:3}] = CalcT1byFAfw2(FFA2DData(:,BestGroup)',ONFAs(BestGroup),TRs(BestGroup));
        
        % Save the first M0 calculated
        First_M0_Map_NFA = RelaxNFA{2};
    else
        % User wanted single M0 for all additional DESPOT calculations
        % Make sure there is only 1 angle in that despot
        if ( numel(FAs) == 1 )
            [RelaxNFA{1:3}] = CalcT1byFAfw2_Single_Angle(FFA2DData(:,BestGroup)',ONFAs(BestGroup),TRs(BestGroup),First_M0_Map_NFA);
        else
            error('-E- There is more than 1 angle in DESPOT when trying to calculate T1');
        end
    end
    
    FRelax4DN=Reshape2DCto4D(RelaxNFA,FBrainMaskS>-1);
    
    if(ExtractFAs)
        NotInBestGroup=setdiff(1:nFAs,BestGroup);
        for ii=1:numel(NotInBestGroup)
            CurFA=NotInBestGroup(ii);
            ONFAs(CurFA)=CalcFAgSigT1M0(FFA2DData(FFBrainMaskSF,CurFA),RelaxNFA{1}(FFBrainMaskSF),RelaxNFA{2}(FFBrainMaskSF),TRsF(CurFA),100,max(0.5,FAsF(CurFA)/2),FAsF(CurFA)*2);
        end
        
        AddToLog(WorkingP,['b_3p' num2str(j)],[num2str(j) ' ONFAsBefore:  ' num2str(ONFAsBefore,'% 8.2f')]);
        AddToLog(WorkingP,['b_3q' num2str(j)],[num2str(j) ' ONFAs:        ' num2str(ONFAs,'% 8.2f')]);
    end
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
            
            
            % Get the set number we are working on
            Set_Num = j;
            
            % Regular DESPOT calculation when this is the first set or
            % when we don't want to use a single M0
            if (Set_Num == 1 || ~Options.Use_Single_M0)
                
                % Calculate T1, PD, RMS for by rest of angles
                [RelaxNewFA{1:3}]=CalcT1byFAfw2(FFA2DDataX(FFBrainMaskSF,BestGroupNewFAs)',NewFAs(BestGroupNewFAs),TRs(BestGroupNewFAs));
                
                
                % Save the first M0 calculated
                First_M0_Map_NewFA = RelaxNewFA{2};
                
            else
                
                % User wanted single M0 for all additional DESPOT calculations
                % Make sure there is only 1 angle in that despot
                if ( numel(FAs) == 1 )
                    [RelaxNewFA{1:3}]=CalcT1byFAfw2_Single_Angle(FFA2DDataX(FFBrainMaskSF,BestGroupNewFAs)',NewFAs(BestGroupNewFAs),TRs(BestGroupNewFAs),First_M0_Map_NewFA);
                else
                    error('-E- There is more than 1 angle in DESPOT when trying to calculate T1');
                end
                
            end
            
            
            
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
        
        if (Options.Calc_Gains_Diff)
            Tmp(:)=SPGRfM(RelaxOFA{1}',RelaxOFA{2}',FAs(i),TRs(i))*T2SFactors(i)*GainFactors(i);
        else
            Tmp(:)=SPGRfM(RelaxOFA{1}',RelaxOFA{2}',FAs(i),TRs(i))*T2SFactors(i);
        end
        
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
    saveas(8181,[WorkingP 'T1andFAestimationSlice' num2str(j) '.png']);
    saveas(8181,[WorkingP 'T1andFAestimationSlice' num2str(j) '.fig']);
    close(8181);
    AddToLog(WorkingP,['yb_5b' num2str(j)],[num2str(j) ' Img T1 fit for Mid slice, Extracted FAs: ' num2str(ONFAs)],['T1andFAestimationSlice' num2str(j) '.png']);
    %%
    Raw2Nii(squeeze(FRelax4DN(:,:,:,1)),[RWorkingP 'T13DNFA.nii'],'float32', MeanFN);
    Raw2Nii(squeeze(FRelax4DN(:,:,:,2)),[RWorkingP 'PD3DNFA.nii'],'float32', MeanFN);
    Raw2Nii(squeeze(FRelax4DN(:,:,:,3)),[RWorkingP 'RMS3DNFA.nii'],'float32', MeanFN);
    
    %% Coregister
    CoregRelaxToMain=true && ~JustT1Mode;
    if(isfield(Options,'CoregRelaxToMain'))
        CoregRelaxToMain=Options.CoregRelaxToMain && ~JustT1Mode;
    end
    if(CoregRelaxToMain && ~strcmp(WhichMean,'Mean 4D'))
        CMeanFN=[WorkingP 'DCEMean.nii'];
        movefile([RWorkingP 'T13DNFA.nii'],[RWorkingP 'BeforeCoreg_T13DNFA.nii']);
        movefile([RWorkingP 'PD3DNFA.nii'],[RWorkingP 'BeforeCoreg_PD3DNFA.nii']);
        movefile([RWorkingP 'RMS3DNFA.nii'],[RWorkingP 'BeforeCoreg_RMS3DNFA.nii']);
        
        T1MatFNs=CoregEstimate([RWorkingP 'BeforeCoreg_T13DNFA.nii'],CMeanFN,CalcForce);
        T1CrgFNs=CoregWrite([RWorkingP 'BeforeCoreg_T13DNFA.nii'],T1MatFNs,CalcForce,Set_Dir,false,CMeanFN);
        movefile(T1CrgFNs,[RWorkingP 'T13DNFA.nii']);
        
        PDMatFNs=CoregEstimate([RWorkingP 'BeforeCoreg_PD3DNFA.nii'],CMeanFN,CalcForce);
        PDCrgFNs=CoregWrite([RWorkingP 'BeforeCoreg_PD3DNFA.nii'],PDMatFNs,CalcForce,Set_Dir,false,CMeanFN);
        movefile(PDCrgFNs,[RWorkingP 'PD3DNFA.nii']);
        
        RMSMatFNs=CoregEstimate([RWorkingP 'BeforeCoreg_RMS3DNFA.nii'],CMeanFN,CalcForce);
        RMSCrgFNs=CoregWrite([RWorkingP 'BeforeCoreg_RMS3DNFA.nii'],RMSMatFNs,CalcForce,Set_Dir,false,CMeanFN);
        movefile(RMSCrgFNs,[RWorkingP 'RMS3DNFA.nii']);
        
        movefile([RWorkingP 'T13DOFA.nii'],[RWorkingP 'BeforeCoreg_T13DOFA.nii']);
        movefile([RWorkingP 'PD3DOFA.nii'],[RWorkingP 'BeforeCoreg_PD3DOFA.nii']);
        movefile([RWorkingP 'RMS3DOFA.nii'],[RWorkingP 'BeforeCoreg_RMS3DOFA.nii']);
        
        T1MatFNs=CoregEstimate([RWorkingP 'BeforeCoreg_T13DOFA.nii'],CMeanFN,CalcForce);
        T1CrgFNs=CoregWrite([RWorkingP 'BeforeCoreg_T13DOFA.nii'],T1MatFNs,CalcForce,Set_Dir,false,CMeanFN);
        movefile(T1CrgFNs,[RWorkingP 'T13DOFA.nii']);
        
        PDMatFNs=CoregEstimate([RWorkingP 'BeforeCoreg_PD3DOFA.nii'],CMeanFN,CalcForce);
        PDCrgFNs=CoregWrite([RWorkingP 'BeforeCoreg_PD3DOFA.nii'],PDMatFNs,CalcForce,Set_Dir,false,CMeanFN);
        movefile(PDCrgFNs,[RWorkingP 'PD3DOFA.nii']);
        
        RMSMatFNs=CoregEstimate([RWorkingP 'BeforeCoreg_RMS3DOFA.nii'],CMeanFN,CalcForce);
        RMSCrgFNs=CoregWrite([RWorkingP 'BeforeCoreg_RMS3DOFA.nii'],RMSMatFNs,CalcForce,Set_Dir,false,CMeanFN);
        movefile(RMSCrgFNs,[RWorkingP 'RMS3DOFA.nii']);
    end
    
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
    
    if ( log(MaxRMSx) > 0)
        subplot(2,2,3);imagesc(mritransform(log(OFARMS)),[0 log(MaxRMSx)]);title('log RMS nominal FAs');
        set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
        subplot(2,2,4);imagesc(mritransform(log(NFARMS)),[0 log(MaxRMSx)]);title('log RMS extracted FAs');
        set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
    else
        subplot(2,2,3);imagesc(mritransform(log(OFARMS)));title('log RMS nominal FAs');
        set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
        subplot(2,2,4);imagesc(mritransform(log(NFARMS)));title('log RMS extracted FAs');
        set(gca,'XTickLabel',[]);set(gca,'YTickLabel',[]);
    end
    
    
    xlabel(num2str(log(MaxRMSx)));
    saveas(91919191,[WorkingP 'T1andFAestimation' num2str(j) '.png']);
    saveas(91919191,[WorkingP 'T1andFAestimation' num2str(j) '.fig']);
    close(91919191);
    
    
    if(ExtractFAs)
        save(MatFN,'ONFAs','TRsF','FAsF','ONFAsBefore');
    else
        save(MatFN,'ONFAs','TRsF','FAsF');
    end
    
    EEBrainMask=imerode(EBrainMask,se);
    EEEBrainMask=imerode(EEBrainMask,se);
    AddToLog(WorkingP,['yb_5c' num2str(j)],[num2str(j) ' Img T1 fit RMS'],['T1andFAestimation' num2str(j) '.png']);
    %% T1 seg removed
    if(false)
        ForSeg=loadniidata(T1Res{CurCGroup,2,1});
        ForSeg(~BrainMask)=NaN;
        ForSegFN=[RWorkingP 'ForSeg.nii'];
        Raw2Nii(ForSeg,ForSegFN,'float32', MeanFN);
        
        T1NFA_SegP=SPM_Segment(ForSegFN,CalcForce,[],false,FMaskFN);
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
        if(ExtractFAs)
            xlabel(num2str([FAs; ONFAs],'% 2.2f'));
        else
            xlabel(num2str([FAs],'% 2.2f'));
        end
        
        saveas(9899,[WorkingP 'T1Seg' num2str(j) '.png']);
        saveas(9899,[WorkingP 'T1Seg' num2str(j) '.fig']);
        close(9899);
        AddToLog(WorkingP,['yb_6' num2str(j)],[num2str(j) ' Img segmentation. Red - GM, Green - WM, Blue - CSF, Magenta - WM for reference, Yellow - CSF for reference.'],['T1Seg' num2str(j) '.png']);
        %%
        Raw2Nii(CSFMask,[WorkingP 'RefAuto' num2str(j) '_CSF_2430.nii'],'float32', MeanFN);
        Raw2Nii(WMMask,[WorkingP 'RefAuto' num2str(j) '_WM_830.nii'],'float32', MeanFN);
    end
end

a=whos;
BigMem={a([a.bytes]>10000).name}';
clear(BigMem{:});
AddToLog(WorkingP,'b_9','Relaxometry finished');
save(RelaxFN);
disp(['DCET1_RelaxForSubjectf finished for ' WorkingP]);
if(~JustT1Mode)
    return;
end
%% Step 6.5 SPM masking
BrainMskFN=[WorkingP 'BrainMask.nii'];
[tmp, MaxFAI]=max(FAsF);

BaselineNoBadSliFN=OtherFACrgFNs{MaxFAI}; % [WorkingP 'BaselineNoBadSli.nii'];
BaselineNoBadSli=loadniidata(BaselineNoBadSliFN);

MeanVol=BaselineNoBadSli;

BaselineNoBadSli(:,:,BadSlicesF2)=NaN;
Raw2Nii(BaselineNoBadSli,BaselineNoBadSliFN,'float32', MeanFN);
% 
% DCEMeanSegP=SPM_Segment(BaselineNoBadSliFN,Force,[],false);
DCEMeanSegQB1=SPM_SegmentWithB1(BaselineNoBadSliFN,CalcForce);

DCEMeanSegP=DCEMeanSegQB1;
C1=loadniidata([DCEMeanSegP 'c1ForSeg.nii'])/256;
C2=loadniidata([DCEMeanSegP 'c2ForSeg.nii'])/256;
C3=loadniidata([DCEMeanSegP 'c3ForSeg.nii'])/256;
if(Options.Mask_thresh>0)
    MinSPMBrainValue=Options.Mask_thresh;
else
    MinSPMBrainValue=0.5;    
end
% % BrainMask=(C1+C2)>MinSPMBrainValue;
% % 
BrainMaskA=(C1+C2+C3)>MinSPMBrainValue;
BrainMaskA=bwfillHoles3Dby2D(BrainMaskA);
se=strel('disk',2,8);
BrainMaskA=imerode(BrainMaskA,se);
BrainMaskA=bwfillHoles3Dby2D(BrainMaskA);
if((Options.Mask_thresh>0) && ~exist(ManMaskFN,'file'))
    BrainMask=BrainMaskA;
    FBrainMask=BrainMask;
end
Raw2Nii(BrainMask,BrainMskFN,'float32',MeanFN);
Raw2Nii(BrainMask,FMaskFN,'float32',MeanFN);

disp('SPM segment finished');
BaselineNoBadSliFN2=[WorkingP 'BaselineNoBadSli2.nii'];

% In unix, run the system cp command with no "-p" because it gives an
% error when the destination is in another computer so source and dest
% files have different owner
if (filesep == '/') % Unix
    system(['cp -f ' BaselineNoBadSliFN ' ' BaselineNoBadSliFN2]);
else  % Windows
    copyfile(BaselineNoBadSliFN,BaselineNoBadSliFN2,'f');    
end

% DCEMeanSegP2=SPM_Segment(BaselineNoBadSliFN2,Force,[],FMaskFN);
DCEMeanSegP2=DCEMeanSegQB1;

BaseSeg3D(:,:,:,1)=loadniidata([DCEMeanSegP2 'c1ForSeg.nii'])/256;
BaseSeg3D(:,:,:,2)=loadniidata([DCEMeanSegP2 'c2ForSeg.nii'])/256;
BaseSeg3D(:,:,:,3)=loadniidata([DCEMeanSegP2 'c3ForSeg.nii'])/256;
BaseSeg3D(:,:,:,4)=loadniidata(FMaskFN);
[Tmp, BaseSeg3DAll]=max(BaseSeg3D(:,:,:,1:3),[],4);
BaseSeg3DAll(~BaseSeg3D(:,:,:,4))=0;
% BaseCleaned=loadniidata([DCEMeanSegP2 'mForSeg.nii']);
BaseCleaned=MeanVol;

se=strel('disk',4,8);
EBrainMask=imerode(FBrainMask,se);
EEBrainMask=imerode(EBrainMask,se);
EEEBrainMask=imerode(EEBrainMask,se);

CSFMask=BaseSeg3DAll==3 & BaseSeg3D(:,:,:,3)>Options.ThreshForRefMasks & EEEBrainMask==1;
WMMask=BaseSeg3DAll==2 & BaseSeg3D(:,:,:,2)>Options.ThreshForRefMasks & EEEBrainMask==1;

if(sumn(WMMask)<100)
    WMMask=BaseSeg3DAll==2 & BaseSeg3D(:,:,:,2)>Options.ThreshForRefMasks*0.9 & EEEBrainMask==1 & ~TooEnhancedForNAWM;
end
BaseSeg3DAllx=BaseSeg3DAll;
BaseSeg3DAllx(CSFMask)=4;
BaseSeg3DAllx(WMMask)=5;

Tmp=max(BrainMask,[],3);
F=find(max(Tmp,[],2));
GoodRows=F(1):F(end);
F=find(max(Tmp,[],1));
GoodCols=F(1):F(end);
GoodSlices=setdiff(1:SDCE(3),BadSlicesF2);

MaxV=median(BaseCleaned(isfinite(BaseCleaned) & BaseCleaned>100))*2;
[q MaxV]=FindDR(BaseCleaned(BrainMask));
for i=1:numel(GoodSlices)
    CurSli=GoodSlices(i);
    I=squeeze(BaseCleaned(:,:,CurSli));
    %         [Tmp MaxV]=FindDR(I);
    ClrM=[1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1];
    IRGB=repmat(min(1,I/MaxV),[1 1 3]);
    for tt=1:3
        BW2 = bwmorph(squeeze(BaseSeg3DAll(:,:,CurSli)==tt),'remove');
        for kk=1:3
            TmpI=IRGB(:,:,kk);
            TmpI(BW2)=ClrM(tt,kk);
            IRGB(:,:,kk)=TmpI;
        end
    end
    for tt=4:5
        BW2 = bwmorph(squeeze(BaseSeg3DAllx(:,:,CurSli)==tt),'remove');
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
saveas(9899,[WorkingP 'BaseSeg'  '.png']);
saveas(9899,[WorkingP 'BaseSeg'  '.fig']);
%%
close(9899);
AddToLog(WorkingP,'ya_2d',['Img segmentation. Red - GM, Green - WM, Blue - CSF, Magenta - WM for reference, Yellow - CSF for reference.'],['BaseSeg'  '.png']);
%%
Raw2Nii(CSFMask,[WorkingP 'RefAuto_Base' '_CSF_2430.nii'],'float32', MeanFN);
Raw2Nii(WMMask,[WorkingP 'RefAuto_Base' '_WM_830.nii'],'float32', MeanFN);

%% B1
T1CleanedFN=[RWorkingP 'Seg_ForSeg'  filesep 'mForSeg.nii'];
% T1UnCleanedFN=[RWorkingP 'Seg_ForSeg'  filesep 'ForSeg.nii'];
T1UnCleanedFN=[RWorkingP 'T13DNFA.nii'];
% RelaxFN=T1CleanedFN;
RelaxFN=T1UnCleanedFN;
% Load .nifti data of the T1 map
CT1=loadniidata(RelaxFN);

% Find B1 by WM ref
D = dir([WorkingP 'RefT1_*.nii']);
if(isempty(D))
    D=dir([WorkingP 'RefAuto_Base_WM*.nii']);
end

Baseline=MeanVol;

RefFN=[WorkingP D(1).name];
AddToLog(WorkingP,'c_20Norm',['T1 Normalization ref file: ' strrep(D(1).name,'_','-')]);
RefTrgVal=str2num(RefFN(find(RefFN=='_',1,'last')+1:end-4));
RefVol=loadniidata(RefFN)>0;
D2 = dir([WorkingP 'WMExMask.nii']);
if(~isempty(D2))
    AddToLog(WorkingP,'c_20NormEx','Using WM exclusion mask!');
    ExMask=loadniidata([WorkingP 'WMExMask.nii'])>0;
    RefVol=RefVol & ~ExMask;
end
RefVolX=RefVol & CT1>100;
% If removing the bad slices, leave us with no remaining voxels, keep them
TmpRefVol = RefVolX;
TmpRefVol(:,:,BadSlicesF2) = false;
if sum(sum(sum(TmpRefVol))) ~= 0
    RefVolX(:,:,BadSlicesF2) = false;
end  
clear TmpRefVol;

%%
% ASK GILAD - what is the purpose of the following?
if(Options.UseN3OnT1)
    %     PseudoB1=sqrt(double(loadniidata(T1Res{1,DoDeviations+1,2}))./double(loadniidata(T1Res{1,DoDeviations+1,1})));
    %     UncleanedT1=loadniidata(T1Res{1,1+Options.ExtractFAs,1});
    UncleanedT1=loadniidata(T1UnCleanedFN);
%     PseudoB1=sqrt(CT1./UncleanedT1);
    ShowSurfFit=false;
    NonLinTransformPower=1;
    NonLinTransform=@(x) log(x.^(1/NonLinTransformPower));
    rNonLinTransform=@(x) exp(x).^NonLinTransformPower;

    FindB1gWM_LS;
    %%
%     if(ShowSurfFit)
    gfig(5001);clf;
    for i=1:numel(GoodSlices)
            CurSli=GoodSlices(i);
%         CurSli=11;
        F=find(I{3}==CurSli);
        FF=find(IF{3}==CurSli);
        [X Y]=meshgrid(1:size(RefVol,2),1:size(RefVol,1));
        X(:)=IF{1}(FF);
        Y(:)=IF{2}(FF);
        gsubplot(numel(GoodSlices),i);
        CurFld=squeeze(NewField(:,:,CurSli));
        surf(X,Y,CurFld,'EdgeColor','none');
        colormap hsv
        alpha(.7)
        hold on;
        if(isempty(F))
%             plot3(I{1}(F), I{2}(F),rNonLinTransform(SolVec(F))/RefTrgVal,'ko');
%             ax=axis;
%             axis([min(I{1}(F))-10 max(I{1}(F))+10 min(I{2}(F))-10 max(I{1}(F))+10 ax(5:6)]);
        else
            plot3(I{1}(F), I{2}(F),rNonLinTransform(SolVec(F))/RefTrgVal,'ko');
            ax=axis;
            axis([min(I{1}(F))-10 max(I{1}(F))+10 min(I{2}(F))-10 max(I{2}(F))+10 ax(5:6)]);
        end
        title(CurSli);
        zlabel([mean(CurFld(:)) std(CurFld(:))]);
    end
    %%
    MaximizeSaveCloseAndAddToLog(WorkingP,'B1gWM','yc_30B1gWMrRMS','$B_1$ given WM.');
    %%
    B1FN=[RWorkingP 'N3B1.nii'];
    Raw2Nii(PseudoB1,B1FN,'float32', MeanFN);
end
%% Setting threshold values

% ASK GILAD - how did he determine those values?
% ASNWER - He uses those values to see if the T1 calculation failed.
%                    He checks whether 1000 voxels or more has a value greater than T1 = 4000.
%                    He used T1=1200 before for the iterations for calculating the REAL T1_base to converge.
%                    T1=40,000 before was just a really high limit (so we won't realty limit the possible results).

% Maximum value of T1 in image
T1Thresh=10000;

% Maximum number of values exceeding the T1 threshold
NaTThresh=1000;

Msk=FBrainMask;

% Go over each slice
for i=1:size(CT1,3)
    CurMsk=squeeze(Msk(:,:,i));
    CurCT1=squeeze(CT1(:,:,i));
    nAT(i)=sumn(CurCT1(CurMsk)>T1Thresh);
end
BadSlicesAgain=find(nAT>NaTThresh);
% If the number of pixels in slice, exceeding T1Thresh is bigger than NaTThresh, mark the slice as NaN.
CT1(:,:,BadSlicesAgain)=NaN;
%%
% ASK GILAD - I don't understand the B1 calculations
% ASNWER       - SHOULD READ THE ARTICLE AND CHECK THAT.
if(Options.UseN3OnT1)
    B1=loadniidata(B1FN);
else
    B1=CT1*0+1;
end

S=sort(CT1(RefVolX));
WhichPercent=0.5;
RefCurVal=S(max(1,floor(numel(S)*WhichPercent)));
% RefCurVal=median(CT1(RefVolX));
kCoeff=RefTrgVal/RefCurVal;
AddToLog(WorkingP,'c_1kc',['kCoeff:' num2str(kCoeff)]);
FinalT1=CT1*kCoeff;
FinalT1FromDESPOT1=FinalT1;


if(~Options.UseN3OnT1)
    kCoeff=1;
end
% kCoeff=1;
B1Nok=B1;
B1=B1*sqrt(kCoeff);
FinalT1FN=[WorkingP 'FinalT1.nii'];
Raw2Nii(FinalT1,FinalT1FN,'float32',MeanFN);

WMD=dir([WorkingP 'RefT1_*WM*.nii']);
if(isempty(WMD))
    WMD=dir([WorkingP 'Ref*.nii']);
    WMD=WMD(strhas({WMD.name},'WM'));
end
RefWMFN=[WorkingP WMD(1).name];
RefWMVol=loadniidata(RefWMFN)>0;
RefWMVal=median(FinalT1(RefWMVol));
AddToLog(WorkingP,'c_2wm',['RefWMVal:' num2str(RefWMVal)],[],2);