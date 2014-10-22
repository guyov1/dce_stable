% script of simulation - to test deconvolution methods
clc;
close all;
clear all;

%---------- R(t):  (high resolution)  --------------
MTT=4; %[sec]
CBF=60; %[ml/100 g/min]
deltaT=0.005; %[sec]
deltaT_DSC=2;
N_time_points=120/deltaT; % The total time will be 2 minutes.
t_vec=0:deltaT:(N_time_points-1)*deltaT;
R_exp_true=exp(-t_vec/MTT);
R_lor_true=1./(1+(pi*t_vec/(2*MTT)).^2);

% figure;plot(R_exp_true,'k');hold on;plot(R_lor_true,'r-.')
% title('examples for R(t)');
% legend('R_e_x_p(t) true','R_l_o_r(t) true');
%----------------------------------------------------



%----------- AIF(t):  (high resolution) -------------
%from Guy/Gilad:
% AIF average population parameters
A1    = 0.809;
A2    = 0.330;

% Changed to shift the boluses to the right
T1    = 0.17046;
T2    = 0.365;
% T1    = 0.27046;
% T2    = 0.465;
% T1    = 0.32046;
% T2    = 0.515;


sig1  = 0.0563;
sig2  = 0.132;
alpha = 1.050;
beta  = 0.1685;
s     = 38.078;
tau   = 0.483;


AIF=AIF_Parker(t_vec/60,A1,sig1,T1,A2,sig2,T2,alpha,beta,s,tau);

% figure;plot(t_vec,AIF,'k');
% title('AIF(t), synthesized by AIF_Parker function');

%--------------------------------------------------


%---------  convolution  AIF(t)*R(t): -------------
fac=1;
% using "filter":
% Ct_signal_Rexp=filter(R_exp_true,sum(R_exp_true),fac*AIF);
% Ct_signal_Rlor=filter(R_lor_true,sum(R_lor_true),fac*AIF);

Ct_signal_Rexp=filter(R_exp_true,1/deltaT,fac*AIF);
Ct_signal_Rlor=filter(R_lor_true,1/deltaT,fac*AIF);

% using "conv":
% Ct_signal_Rexp=conv(R_exp_true,fac*AIF);
% Ct_signal_Rlor=conv(R_lor_true,fac*AIF);


% figure;plot(Ct_signal_Rexp,'k');hold on;plot(Ct_signal_Rlor,'r-.')
% title('C(t) = CBF(AIF(t)*R(t))');
% legend('C(t) of voxel with R_e_x_p(t)','C(t) of voxel with R_l_o_r(t)');
%---------------------------------------------------



%---------- Sampling like in real MRI data: -------------
delta_sampling_time_points=deltaT_DSC/deltaT;
Ct_signal_samp_Rexp=downsample(Ct_signal_Rexp,delta_sampling_time_points);
Ct_signal_samp_Rlor=downsample(Ct_signal_Rlor,delta_sampling_time_points);
AIF_samp=downsample(AIF,delta_sampling_time_points);
R_exp_true_samp=downsample(R_exp_true,delta_sampling_time_points);
R_lor_true_samp=downsample(R_lor_true,delta_sampling_time_points);

% figure;plot(Ct_signal_samp_Rexp,'k');hold on;plot(Ct_signal_samp_Rlor,'r-.');
% title(['C(t) after sampled in real time points (every ',num2str(deltaT_DSC),'[sec])']);
%------------------------------------------------------------


%----------  Add gaussian noise  ----------------------------
SNR=20;
signal_power_Rexp=sum(Ct_signal_samp_Rexp.^2)/(length(Ct_signal_samp_Rexp)*deltaT_DSC);
signal_power_Rlor=sum(Ct_signal_samp_Rlor.^2)/(length(Ct_signal_samp_Rlor)*deltaT_DSC);

noise_power_Rexp=signal_power_Rexp*10^(-SNR/10);
noise_power_Rlor=signal_power_Rlor*10^(-SNR/10);

Ct_noise_Rexp=sqrt(noise_power_Rexp).*randn(1,length(Ct_signal_samp_Rexp));
Ct_noise_Rlor=sqrt(noise_power_Rlor).*randn(1,length(Ct_signal_samp_Rlor));

Ct_Rexp=Ct_signal_samp_Rexp+Ct_noise_Rexp;
Ct_Rlor=Ct_signal_samp_Rlor+Ct_noise_Rlor;

figure;plot(Ct_Rexp,'k');hold on;plot(Ct_Rlor,'r-.');
title(['C(t) after noise was added. SNR=',num2str(SNR),'[dB]']);
%------------------------------------------------------------


%----------- convert AIF to matrix:  -------------------------
[AIFmat] = AIF2mat(AIF_samp,deltaT_DSC,'linear');



LineSpecs={'k-','k-.','r-','r-.'};

for R_ind=1:1
    if R_ind==1
        Ct=Ct_Rexp;
        Rt_true_samp=R_exp_true_samp;
    else
        Ct=Ct_Rlor;
        Rt_true_samp=R_lor_true_samp;
    end
    
    %------------------ Solve using normal SVD ----------------
    
    [U,S,V] = svd(AIFmat);
    
    %define threshold (% of the max of S):
    PSVD_perc=15;
    PSVD=max(diag(S))*PSVD_perc/100;
    PSVD_inds=find(diag(S)<PSVD);
    
    % calc W=1/S and cut small values (small in S):
    W_diag=1./diag(S);
    W_diag(PSVD_inds)=0;
    W=diag(W_diag);
    
    %estimate R(t):
    Rt_SVD=V*W*U'*Ct(:);
    % figure:
    figure;plot(Rt_SVD,LineSpecs{2*R_ind-1});
    hold on;
    plot(Rt_true_samp,LineSpecs{2*R_ind});
    title(['reconstructed R(t) vs. true R(t), Using normal SVD method.   PSVD = ',num2str(PSVD_perc),'% of max S']);
    legend('R(t) reconstructed','R(t) ground truth');
    
    %------------------ Solve using circular SVD (oSVD & cSVD)----------------
    L=2*length(AIF_samp); % L=2N, for circular SVD
    % in circular SVD we zero-padd the AIF and Ct, and build a new AIF-mat.
    % The lower triangle is built the same way. The upper triangle is
    % built s.t. the whole matrix is circular
    AIF_zero_padd=zeros(L,1);
    Ct_zero_padd=zeros(L,1);
    AIF_zero_padd(1:length(AIF_samp))=AIF_samp;
    Ct_zero_padd(1:length(Ct))=Ct;
    % build the lower triangle of the new AIF matrix:
    AIFmat_circular=AIF2mat(AIF_zero_padd,deltaT_DSC,'linear');
    % build the upper triangle, according to 2003 paper:
    for ii=1:L-1
        AIFmat_circular(1,ii+1:L)=fliplr(AIFmat_circular(ii+1:L,1)');
    end
    
    [U_circ,S_circ,V_circ] = svd(AIFmat_circular);
    
    %define threshold for cSVD (% of the max of S):
    PSVD_circ_cSVD_perc=20;
    PSVD_circ_cSVD=max(diag(S_circ))*PSVD_circ_perc/100;
    PSVD_circ_cSVD_inds=find(diag(S_circ)<PSVD_circ);
    
    % calc W=1/S and cut small values (small in S):
    W_circ_diag_cSVD=1./diag(S_circ); %for oSVD (WITH regularization)
    W_circ_diag_cSVD(PSVD_circ_cSVD_inds)=0;
    W_circ_cSVD=diag(W_circ_diag_cSVD);
    
    %estimate R(t):
    %cSVD:
    Rt_cSVD=V_circ*W_circ_cSVD*U_circ'*Ct_zero_padd(:);
    
    % figure:
    figure;plot(Rt_cSVD(1:length(AIF_samp)),LineSpecs{2*R_ind-1});
    hold on;
    plot(Rt_true_samp(1:length(AIF_samp)),LineSpecs{2*R_ind});
    title(['reconstructed R(t) vs. true R(t), Using oSVD (circular) method.   PoSVD = ',num2str(PSVD_circ_perc),'% of max S']);
    legend('R(t) reconstructed','R(t) ground truth');
    
    %------------------ Solve using Tikhonov regularization ---------
    
    
    L1=zeros(size(AIFmat,1)-1,size(AIFmat,2));
    zero_vec=zeros(size(Ct));
    
    L1_main_diagonal=(-1)*ones(size(L1,2),1);
    L1_second_diagonal=ones(size(L1,2)-1,1);
    
    L1a=diag(L1_main_diagonal,0);
    L1b=diag(L1_second_diagonal,1);
    
    L1=L1a(1:end-1,:)+L1b(1:end-1,:);
    
    %make L1 a square matrix (by adding a row), for inverting it and use it for standard form
    L1_sq=zeros(size(L1,2));
    L1_sq(1:end-1,:)=L1;
    L1_sq(end,end)=-1;
    
    %make L1 closer to derivative by dividing in the time diff between
    %points:
    L1_sq=L1_sq/deltaT_DSC;
    
    %--- generalized SVD:
    [U,V,X,C,S] = gsvd(AIFmat,L1);
    mu=diag(S);
    s=diag(C);
    sigma=gsvd(AIFmat,L1);
    % sigma=sqrt(diag(C'*C)./diag(S'*S));
    %     sm=[reshape(sigma(1:end-1),length(sigma)-1,1) mu(:)];
    sm=[s(1:end-1) mu(:)];
    
    % plot the L-curve:
    %     [reg_corner,rho,eta,reg_param]=l_curve(U,s,Ct(:),'Tikh',L1,V);
    %     [reg_corner,rho,eta,reg_param]=l_curve(U,sm,Ct(:),'Tikh');
    % Solve R(t) by Tikhonov regulrization:
    %     [Rt,rho,eta] = tikhonov(U,sm,X,Ct(:),sqrt(reg_corner));
    
    %-----------------------
    
    %--- Tranforming to standard form (after squaring L1):
    AIFmat_sf=AIFmat/(L1_sq);
    
    [U,W,V] = svd(AIFmat_sf);
    s=diag(W);
    
    %     [reg_corner,rho,eta,reg_param]=l_curve(U,s,Ct(:),'Tikh');
    %     [Rt_sf,rho,eta] = tikhonov(U,s,V,Ct(:),sqrt(reg_corner));
    %     Rt=(L1_sq)\Rt_sf;
    
    % finding the best lambda and best total norm:
    lambda_vec=0.1:0.1:100;
    lambda_best=0;
    norm_best=inf;
    for lambda_ind=1:length(lambda_vec)
        lambda=lambda_vec(lambda_ind);
        [Rt_sf,rho,eta] = tikhonov(U,s,V,Ct(:),lambda);
        Rt=(L1_sq)\Rt_sf;
        %calc norm without standart form
        norms_Ax_b(lambda_ind)=norm(AIFmat*Rt-Ct(:))^2;
        norms_Lx(lambda_ind)=norm(L1*Rt)^2;;
        norm_Ax_b_Lx=norms_Ax_b(lambda_ind)+lambda^2*norms_Lx(lambda_ind);
        norms_Ax_b_Lx(lambda_ind)=norm_Ax_b_Lx;
        if norm_Ax_b_Lx<norm_best
            norm_best=norm_Ax_b_Lx;
            lambda_best=lambda;
        end
    end
    
    %plot L curve:
    figure;
    loglog(sqrt(norms_Ax_b),sqrt(norms_Lx));
    
    lambda_defacto=0.1;%(reg_corner);
    [Rt_sf,rho,eta] = tikhonov(U,s,V,Ct(:),lambda_defacto);
    Rt=(L1_sq)\Rt_sf;
    data_norm=norm(AIFmat*Rt-Ct(:))^2;
    reg_norm=norm(L1*Rt)^2;
    
    figure;plot(Rt,LineSpecs{2*R_ind-1});
    hold on;
    plot(Rt_true_samp,LineSpecs{2*R_ind});
    title(['reconstructed R(t) vs. true R(t).    \lambda=',num2str(lambda_defacto),'.  ||Ax-b||^2=',num2str(data_norm),'. ||Lx||^2=',num2str(reg_norm),'.  ||Ax-b||^2+\lambda^2||Lx||^2=',num2str(data_norm+lambda_defacto^2*reg_norm),'.   ||Rt_r_e_c-Rt_t_r_u_e||^2 = ',num2str(norm(Rt-Rt_true_samp(:))^2)]);
    legend('R(t) reconstructed','R(t) ground truth');
    
    a=5;
    
end
%-------------------------------------------------------------


%---------- deconvolution  (no noise added to AIF so far) -----------


%-------------------------------------------------------------