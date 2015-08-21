function [retunStruct] = chooseAifForRealData(Sim_Struct, CTC2D, Art_Mask, Vein_Mask, Msk2, Output_directory, AIFFindData_mat)
%chooseAifForRealData Summary of this function goes here

time_vec_minutes = Sim_Struct.time_vec_minutes;
Correct_PVE      = Sim_Struct.Correct_PVE;

if exist('AIFFindData_mat','var')
    load(AIFFindData_mat);
    
    % AIF_Parker8t=@(x,t) AIF_Parkerg2( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8))*x(2);
    AIF_Parker9t    = @(x,t)AIF_Parkerg3( t,1,x(3),x(1),x(5),x(6),x(1)+x(4),x(7),x(8),x(9),max(time_vec_minutes))*x(2);
    % Create Parker's AIF
    AIF_parametric = AIF_Parker9t(OutAIFParam,time_vec_minutes);
end

% Artreries and Veins masks
ICA_Art3D   = loadniidata(Art_Mask);
ICA_Art2D   = Reshape4d22d(ICA_Art3D, Msk2);
% Indices in 2D for Arteris/veins
ICA_Art2D_indices  = find(ICA_Art2D>0);

if exist('Vein_Mask','var')
    ICA_Vein3D         = loadniidata(Vein_Mask);
    ICA_Vein2D         = Reshape4d22d(ICA_Vein3D,Msk2);
    ICA_Vein2D_indices = find(ICA_Vein2D>0);
end

AIF_estimated_ICA  = transpose(mean(transpose(CTC2D(ICA_Art2D_indices,:)),2));

if exist('Vein_Mask','var')
    Vein_estimated_ICA = transpose(mean(transpose(CTC2D(ICA_Vein2D_indices,:)),2));
    % Coorect AIF scale to avoid Partial Volume Effect
    if Correct_PVE
        [AIF_estimated_ICA_PVE_correct, Scale_Factor] = CorrectPVE(AIF_estimated_ICA, Vein_estimated_ICA, Sim_Struct);
    end
end

%% Plot Arteris/veins Ct's and their averages
fig_num            = figure;

subplot(1,2,1);
if exist('AIFFindData_mat','var')
    plot(time_vec_minutes,AIF_parametric,time_vec_minutes,AIF_parametric,'*');
end
title('Functional AIF');
xlabel('Time [Min]');
ylabel('Amplitude');
subplot(1,2,2);
hold on;
plot(time_vec_minutes,transpose(CTC2D(ICA_Art2D_indices,:)),'r--');
h1 = plot(time_vec_minutes,mean(transpose(CTC2D(ICA_Art2D_indices,:)),2),'r','LineWidth',4);

if exist('Vein_Mask','var')
    plot(time_vec_minutes,transpose(CTC2D(ICA_Vein2D_indices,:)),'b--');
    h2 = plot(time_vec_minutes,mean(transpose(CTC2D(ICA_Vein2D_indices,:)),2),'b','LineWidth',4);
    if Correct_PVE
        h3 = plot(time_vec_minutes,AIF_estimated_ICA_PVE_correct,'k','LineWidth',4);
    end
end

title('Arteris/Veins Cts and their averages');
xlabel('Time [Min]');
ylabel('Amplitude');
if exist('Vein_Mask','var')
    if Correct_PVE
        legend([h1 h2 h3],'Mean Artery','Mean Vein', 'Mean Artery - Corrected PVE');
    else
        legend([h1 h2],'Mean Artery','Mean Vein');
    end
else
    legend([h1],'Mean Artery');
end

hold off;

%% Add to Log
LogFN = [Output_directory 'Log.mat'];
if exist(LogFN,'file')
    delete(LogFN)
end
SN = 'Report';
Log.idx_000={['\\title{' SN '}\r\n\\maketitle\r\n']};
save(LogFN,'Log');
gprint(fig_num,[Output_directory 'InputAIFs.png']);
%gprint(fig_num,'Run_Output/InputAIFs.png');

AddToLog(Output_directory,'idx_001','\\subsection*{Possible AIF Inputs}');
AddToLog(Output_directory,'idx_002','InputAIFs','InputAIFs.png');
%% Assign to output
retunStruct                        = struct();
if Correct_PVE
    retunStruct.AIF_estimated_ICA     = AIF_estimated_ICA_PVE_correct;    
else
    retunStruct.AIF_estimated_ICA     = AIF_estimated_ICA;
end

if exist('AIF_parametric','var')
    retunStruct.AIF_parametric    = AIF_parametric;
end

end

