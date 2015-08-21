function [Rt] = SVD_solve_voxel(AIF,Ct,deltaT,Rt_type,deconv_methods)

% This function solve (Ab=c) using SVD for b, (as shown in the paper of
% Ostergard 1996), for a SINGLE voxel. 
%   A is matrix built from AIF, 
%   Ct=c=concentration curve at a voxel, Rt_type='constant' or 'linear' (assumption on the
%       difference between 2 time points),
%   b=Rt (the desired solution).
%
sSVD_en=deconv_methods.sSVD.en;
sSVD_th=deconv_methods.sSVD.th;
cSVD_en=deconv_methods.cSVD.en;
cSVD_th=deconv_methods.cSVD.th;
oSVD_en=deconv_methods.oSVD.en;
oSVD_OI=deconv_methods.oSVD.OI;
tikh_en=deconv_methods.tikhonov.en;

AIF=AIF(:);
Ct=Ct(:);

Rt_sSVD=zeros(size(Ct));
Rt_cSVD=zeros(size(Ct));
Rt_oSVD=zeros(size(Ct));
Rt_tikhonov=zeros(size(Ct));

if length(AIF)~=length(Ct)
    error('AIF and concentration curve are not the same length');
end

%------ standard SVD -----------------------%
if sSVD_en 
    AIFmat=AIF2mat(AIF,deltaT,Rt_type);
    [U,S,V] = svd(AIFmat); % A=U*W*V'
    PSVD=max(diag(S))*sSVD_th/100;
    PSVD_inds=find(diag(S)<PSVD);
    % calc W=1/S and cut small values (small in S):
    W_diag=1./diag(S);
    W_diag(PSVD_inds)=0;
    W=diag(W_diag);
    Rt.sSVD=V*W*U'*Ct(:);
end


%------ circular SVD -----------------------%
if cSVD_en || oSVD_en  % Circular SVD
    % in circular SVD we zero-padd the AIF and Ct, and build a new AIF-mat.
    % The lower triangle is built the same way. The upper triangle is
    % built s.t. the whole matrix is circular
    L=2*length(AIF); % L=2N, for circular SVD
    AIF_zero_padd=zeros(L,1);
    Ct_zero_padd=zeros(L,1);
    AIF_zero_padd(1:length(AIF))=AIF;
    Ct_zero_padd(1:length(Ct))=Ct;
    % build the lower triangle of the new AIF matrix:
    AIFmat_circular=AIF2mat(AIF_zero_padd,deltaT,Rt_type);
    % build the upper triangle, according to 2003 paper:
    for ii=1:size(AIFmat_circular,2)-1
        val_to_add_in_diagonal=AIFmat_circular(size(AIFmat_circular,1)-ii+1,1);
        AIFmat_circular=AIFmat_circular+diag(val_to_add_in_diagonal*ones(size(AIFmat_circular,1)-ii,1),ii);
    end
    AIFmat=AIFmat_circular;
    %Now continue in the usual svd flow:
    [U,S,V] = svd(AIFmat_circular); % A=U*W*V'
    % in cSVD use the TH provided by the user.
    % We also do it in oSVD for initial guess of R(t)
    if cSVD_en || oSVD_en
        PSVD=max(diag(S))*cSVD_th/100;
        PSVD_inds=find(diag(S)<PSVD);
        % calc W=1/S and cut small values (small in S):
        W_diag=1./diag(S);
        W_diag(PSVD_inds)=0;
        W=diag(W_diag);
        Rt.cSVD=V*W*U'*Ct_zero_padd(:);
    end
    % in oSVD find the best TH, using Oscillation Index (OI) provided by
    % user:
    if oSVD_en
        Rt_vec=Rt.cSVD(1:length(Ct));
        Rt_OI=sum(abs(Rt_vec(3:end)-2*Rt_vec(2:end-1)+Rt_vec(1:end-2)))/max(Rt_vec)/length(Ct);
        while Rt_OI>oSVD_OI
            cSVD_th=cSVD_th*1.2;
            PSVD=max(diag(S))*cSVD_th/100;
            PSVD_inds=find(diag(S)<PSVD);
            % calc W=1/S and cut small values (small in S):
            W_diag=1./diag(S);
            W_diag(PSVD_inds)=0;
            W=diag(W_diag);
            Rt.cSVD=V*W*U'*Ct_zero_padd(:);
            Rt_vec=Rt.cSVD(1:length(Ct));
            Rt_OI=sum(abs(Rt_vec(3:end)-2*Rt_vec(2:end-1)+Rt_vec(1:end-2)))/max(Rt_vec)/length(Ct);
        end
        Rt.oSVD=Rt_vec;
    end
end

%------ Tikhonov regularization -----------------------%

if tikh_en
    if ~exist('AIFmat')
        AIFmat=AIF2mat(AIF,deltaT,Rt_type);
    end
    % build L1 matrix:
    N=length(Ct);
    P=N-1;
    L1_mat=zeros(P,N);
    L1_mat(:,1:end-1)=-1*eye(P);
    L1_mat(:,2:end)=L1_mat(:,2:end)+eye(P);
    % transformation to standard form, for direct method
%     [K R]=qr(L1_mat.');
%     Kp=K(:,1:end-1);
%     Rp=R(1:end-1,:);
%     Ko=K(:,end);
%     [H T]=qr(AIFmat*Ko);
%     Ho=H(:,1);
%     Hq=H(:,2:end);
%     To=T(1,:);
%     AIFmat_standard=Hq.'*AIFmat*L1_mat';
%     Ct_standard=Hq.'*Ct;
    
    [AIFmat_s,Ct_s,L1_mat_s,K,M] = std_form(AIFmat,L1_mat,Ct);
    [U,S,V] = svd(AIFmat_s);
    %%%%%%%% regularization part, elimintaing small values in S, if needed:
    
    %%%%%%%%%%
    % Solve using SVD and then Tikhonov:
    s=diag(S);
    [reg_corner,rho,eta,reg_param] = l_curve(U,s,Ct_s,'Tikh');
    
    % finding the best lambda and best total norm:
    lambda_vec=0.0001:0.0001:1;
    lambda_best=0;
    norm_best=inf;
    for lambda_ind=1:length(lambda_vec)
        lambda=lambda_vec(lambda_ind);
        [Rt_sf,rho,eta] = tikhonov(U,s,V,Ct_s(:),lambda);
        norms_Ax_b_s(lambda_ind)=norm(AIFmat_s*Rt_sf-Ct_s(:))^2;
        norms_x_s(lambda_ind)=norm(Rt_sf)^2;
%         Rt=gen_form(L1_mat_s,Rt_sf,AIFmat,Ct,K,M);%(L1_sq)\Rt_sf;
        %calc norm without standart form
%         norms_Ax_b(lambda_ind)=norm(AIFmat*Rt-Ct(:))^2;
%         norms_Lx(lambda_ind)=norm(L1_mat*Rt)^2;;
        norm_Ax_b_Lx_s=norms_Ax_b_s(lambda_ind)+lambda^2*norms_x_s(lambda_ind);
        norms_Ax_b_Lx_s(lambda_ind)=norm_Ax_b_Lx_s;
        if norm_Ax_b_Lx_s<norm_best
            norm_best=norm_Ax_b_Lx_s;
            lambda_best=lambda;
        end
    end
    
    lambda=sqrt(reg_corner);
    [Rt_lambda_s,rho,eta] = tikhonov(U,s,V,Ct_s,lambda);
    
    %Return to the original R(t) in the general form:
    Rt.Tikh = gen_form(L1_mat_s,Rt_lambda_s,AIFmat,Ct,K,M);
    Rt_OI=sum(abs(Rt.Tikh(3:end)-2*Rt.Tikh(2:end-1)+Rt.Tikh(1:end-2)))/max(Rt.Tikh)/length(Ct);
    a=5;
end


% Rt_voxel = [Rt_sSVD Rt_cSVD Rt_oSVD];


% [U,W,V] = svd(AIFmat); % A=U*W*V'
% W_inv=inv(W);
% W_inv_after_proc=process_W_SVD(W_inv,th_perc);  % here we can "regularize" the problem by eliminating elements from diagonal
% Rt_voxel_circular= V*W_inv_after_proc*U'*Ct_zero_padd(:); %  b= V * 1/W * U' * c

% AIFmat_inv=inv(AIFmat);
% [Va,Wa,Ua]=svd(AIFmat_inv); % A^-1 = V*W*U'
% Wa_after_proc=process_W_SVD(Wa);
% Rt_voxel_a=Va*Wa_after_proc*Ua'*Ct;  % b=V*W*U'*c


% if tikhonov  % with regularization, using Generalized SVD and Tikhonov Regularization
%     AIFmat=AIF2mat(AIF,deltaT,Rt_type);
%     
%     L1=zeros(size(AIFmat,1)-1,size(AIFmat,2));
%     zero_vec=zeros(size(Ct));
%     
%     L1_main_diagonal=(-1)*ones(size(L1,2),1);
%     L1_second_diagonal=ones(size(L1,2)-1,1);
%     
%     L1a=diag(L1_main_diagonal,0);
%     L1b=diag(L1_second_diagonal,1);
%     
%     L1=L1a(1:end-1,:)+L1b(1:end-1,:);
%     
%     % generalized SVD:
%     [U,V,X,C,S] = gsvd(AIFmat,L1);
%     mu=diag(S);
%     s=diag(C);
%     sigma=gsvd(AIFmat,L1);
%     sm=[s(1:end-1) mu(:)];
%  
%      % plot the L-curve:
%     [reg_corner,rho,eta,reg_param]=l_curve(U,sm,Ct(:),'Tikh');
%     % Solve R(t) by Tikhonov regulrization:
%     [Rt_voxel_tikhonov,rho,eta] = tikhonov(U,sm,X,Ct(:),sqrt(reg_corner));
% end


% figure;plot(Rt_voxel(1:70),'g');hold on;plot(Rt_voxel_circular(1:70),'k-.');plot(Rt_voxel_tikhonov(1:70),'r--');
% title(['TH for SVD=',num2str(th_perc),'%'],'FontSize',12);
% legend('SVD','Circular SVD','gSVD and Tikhonov regularization');
