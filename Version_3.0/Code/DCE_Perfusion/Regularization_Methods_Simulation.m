function [ ridge_regression_result, b_spline_result, b_spline_result_1st_deriv, b_spline_result_2nd_deriv, b_PCA_larss_result, b_PCA_result_1st_deriv, b_PCA_result_2nd_deriv, idx_fig ] = ...
    Regularization_Methods_Simulation( Sim_Ct_T, Sim_Ct_T_noise, Conv_Matrix, Conv_Matrix_no_noise, time_vec_minutes, lambda_vec, normalize,...
    min_interval, B_mat, B_PCA, plot_L_curve, idx_fig, filter_type, Derivative_Time_Devision, plot_flag, RealData_Flag )

% Regularization_Methods_Simulation - Implement deconvolution by regularization
%
% Formalizing the problem:
%
%   Without splines:
%
%   A*R_vec = C_vec
%   min { ||A*R_vec - C_vec||^2  + lambda^2||L*R_vec ||^2 }
%
%   where -
%            A     - Conv_Matrix
%            R_vec - Model vector ( exp. for Larson - F*IRF(t) )
%            C_vec - Ct(t) - Sim_Ct_T_noise - Tissue concentration
%            L     - Some matrix ( Identity/1st Derivative/2nd Derivative)
%
%   With splines:
%   Force the solution R_vec to be in a splines form. Hence,
%   R_vec = B*V_vec
%
%   where -
%            B     - Spline matrix
%            V_vec - Polynomial coefficients for B splines
%
%   So the new minimzation problem, takes the form:
%
%   A*B*V_vec = C_vec
%   If we define A*B = D, we can rewrite:
%   D*V_vec = C_vec
%   min { ||D*V_vec - C_vec||^2  + lambda^2||L*V_vec ||^2 }

% Emphasize IRF(t1) from which we estimate the flow
%Conv_Matrix(1,:)  = 10*Conv_Matrix(1,:);
%Sim_Ct_T_noise(1) = 10*Sim_Ct_T_noise(1);

x0     = zeros(size(Sim_Ct_T_noise));
knots0 = zeros(size(B_mat,2),1);

if ~RealData_Flag
    [ ridge_regression_result ] = Regularized_Sol( Conv_Matrix      , Sim_Ct_T_noise, x0    , lambda_vec(1), normalize, 1);
    [ knots_coeff ]             = Regularized_Sol( Conv_Matrix*B_mat, Sim_Ct_T_noise, knots0, lambda_vec(2), normalize, 0);
    
    % Using derivative regularization constraint
    num_points   = length(Sim_Ct_T_noise);
end

% --------------------------------------------------------------------------------
% First derivative
%deriv_matrix = (1/min_interval)*toeplitz([-1,zeros(1,num_points-1)],[[-1 1],zeros(1,num_points-2)]);
if (Derivative_Time_Devision)
    deriv_matrix_1st = (1/min_interval)*toeplitz([-1,zeros(1,size(B_mat,2)-1)],[[-1 1],zeros(1,size(B_mat,2)-2)]);
else
    deriv_matrix_1st = toeplitz([-1,zeros(1,size(B_mat,2)-1)],[[-1 1],zeros(1,size(B_mat,2)-2)]);
end

%deriv_matrix(50:end,:) = 500*deriv_matrix(50:end,:);
% Second derivative
deriv_matrix_2nd =  deriv_matrix_1st * deriv_matrix_1st ;
%deriv_matrix_2 =  (1/min_interval)*toeplitz([1,zeros(1,num_points-1)],[[1 -2 1],zeros(1,num_points-3)]);


% Display Knots and its derivatives when no noise in encountered
if (plot_flag && ~RealData_Flag)
    
    [ knots_coeff_4_display ]       = Regularized_Sol( Conv_Matrix_no_noise*B_mat , Sim_Ct_T, knots0, 0, 1, 0);
    [ knots_coeff_4_display_noisy ] = Regularized_Sol( Conv_Matrix*B_mat          , Sim_Ct_T_noise, knots0, 0, 1, 0);
    
    
    b_spline_result_4_display       = B_mat*knots_coeff_4_display;
    b_spline_result_4_display_noisy = B_mat*knots_coeff_4_display_noisy;
    
    % Keep positive
    b_spline_result_4_display(b_spline_result_4_display<0)             = 0;
    b_spline_result_4_display_noisy(b_spline_result_4_display_noisy<0) = 0;
    
    fig_num = figure;
    subplot(4,1,1);
    hold on;
    plot(time_vec_minutes,Sim_Ct_T,'g*');
    plot(time_vec_minutes,Conv_Matrix_no_noise * b_spline_result_4_display,'r-');
    plot(time_vec_minutes,Conv_Matrix          * b_spline_result_4_display_noisy,'mo');
    hold off;
    title('B-Spline estimation with and without noise');
    subplot(4,1,2);
    hold on;
    plot(knots_coeff_4_display,'go');
    %plot(knots_coeff_4_display_noisy,'r*');
    hold off;
    title('Knots (and noisy one)');
    subplot(4,1,3);
    hold on;
    plot(deriv_matrix_1st*knots_coeff_4_display,'go');
    %plot(deriv_matrix*knots_coeff_4_display_noisy,'r*');
    hold off;
    title('Knots 1st derivative (and noisy one)');
    subplot(4,1,4);
    hold on;
    plot(deriv_matrix_2nd*knots_coeff_4_display,'go');
    %plot(deriv_matrix_2*knots_coeff_4_display_noisy,'r*');
    hold off;
    title('Knots 2nd derivative (and noisy one)');
    % Print result to PDF
    Fig_name = [filter_type '_Knots_curves.png'];
    Title    = [filter_type ' Knots-Curves'];
    subTitle = '';
    [idx_fig] = Print2Pdf(fig_num, idx_fig, Fig_name, './Run_Output/',Title, subTitle);
    
end


%
%A_new_1 = Conv_Matrix*inv(deriv_matrix)  *B_mat;
%A_new_2 = Conv_Matrix*inv(deriv_matrix_2)*B_mat;
if ~RealData_Flag
    A_new_1 = Conv_Matrix*B_mat*inv(deriv_matrix_1st);
end
A_new_2 = Conv_Matrix*B_mat*inv(deriv_matrix_2nd);

if ~RealData_Flag
    PCA_deriv_matrix_1 = deriv_matrix_1st;
    PCA_deriv_matrix_2 = deriv_matrix_2nd;
    %PCA_deriv_matrix = deriv_matrix;
    A_PCA              = Conv_Matrix * B_PCA;
    A_PCA_1            = Conv_Matrix * B_PCA * inv(PCA_deriv_matrix_1);
    A_PCA_2            = Conv_Matrix * B_PCA * inv(PCA_deriv_matrix_2);
end

% --------------------------------------------------------------------------------

if ~RealData_Flag
    
    [ knots_coeff_PCA ]   = Regularized_Sol( A_PCA      , Sim_Ct_T_noise, knots0, lambda_vec(5), normalize, 0 );
    [ knots_coeff_PCA_1 ] = Regularized_Sol( A_PCA_1    , Sim_Ct_T_noise, knots0, lambda_vec(6), normalize, 0 );
    [ knots_coeff_PCA_2 ] = Regularized_Sol( A_PCA_2    , Sim_Ct_T_noise, knots0, lambda_vec(7), normalize, 0 );
    [ knots_coeff_2 ]     = Regularized_Sol( A_new_1    , Sim_Ct_T_noise, knots0, lambda_vec(3), normalize, 0);
    
    inv_matrix_PCA_1st_deriv  = B_PCA * inv(PCA_deriv_matrix_1);
    inv_matrix_PCA_2nd_deriv  = B_PCA * inv(PCA_deriv_matrix_2);
    inv_matrix_1st_deriv      = B_mat * inv(deriv_matrix_1st) ;
    
    b_PCA_larss_result        = B_PCA        * knots_coeff_PCA;
    b_PCA_result_1st_deriv    = inv_matrix_PCA_1st_deriv * knots_coeff_PCA_1;
    b_PCA_result_2nd_deriv    = inv_matrix_PCA_2nd_deriv * knots_coeff_PCA_2;
    b_spline_result           = B_mat        * knots_coeff;
    b_spline_result_1st_deriv = inv_matrix_1st_deriv * knots_coeff_2;

end

[ knots_coeff_3 ]         = Regularized_Sol( A_new_2    , Sim_Ct_T_noise, knots0, lambda_vec(4), normalize, 0 );
inv_matrix_2nd_deriv      = B_mat * inv(deriv_matrix_2nd);
b_spline_result_2nd_deriv = inv_matrix_2nd_deriv * knots_coeff_3;

% --------------------------------------------------------------------------------

% %%% Testing GSVD  (get_l -> derivative matrices)
% %[UU_1,sm_1,XX_1] = cgsvd(Conv_Matrix*B_mat,get_l(28,1));
% %[UU_2,sm_2,XX_2] = cgsvd(Conv_Matrix*B_mat,get_l(28,2));
%
% deriv_matrix_new   = (1/min_interval)*toeplitz([-1,zeros(1,size(B_mat,2)-1)],[[-1 1],zeros(1,size(B_mat,2)-2)]);
% deriv_matrix_new_2 =  deriv_matrix_new*deriv_matrix_new;
%
% [UU_1,sm_1,XX_1] = cgsvd(Conv_Matrix*B_mat,deriv_matrix_new);
% [UU_2,sm_2,XX_2] = cgsvd(Conv_Matrix*B_mat,deriv_matrix_new_2);
%
% %knots_coeff_2 = lsqi(UU_1,sm_1,XX_1,Sim_Ct_T_noise,4);
% %knots_coeff_3 = lsqi(UU_2,sm_2,XX_2,Sim_Ct_T_noise,4);
%
% knots_coeff_2 = tikhonov(UU_1,sm_1,XX_1,Sim_Ct_T_noise,lambda_vec(3));
% knots_coeff_3 = tikhonov(UU_2,sm_2,XX_2,Sim_Ct_T_noise,lambda_vec(4));
%
% b_spline_result_1st_deriv = B_mat*knots_coeff_2;
% b_spline_result_2nd_deriv = B_mat*knots_coeff_3;

% --------------------------------------------------------------------------------

% Remove zeros
if ~RealData_Flag
    b_PCA_larss_result        = max(b_PCA_larss_result,0);
    b_PCA_result_1st_deriv    = max(b_PCA_result_1st_deriv,0);
    b_PCA_result_2nd_deriv    = max(b_PCA_result_2nd_deriv,0);
    b_spline_result           = max(b_spline_result,0);
    b_spline_result_1st_deriv = max(b_spline_result_1st_deriv,0);
end
b_spline_result_2nd_deriv = max(b_spline_result_2nd_deriv,0);

%% Create L curve
if (plot_L_curve )
    
    %lambda_vec_string = '0.000001:0.010:3';
    %lambda_vec_string = '0.3:0.1:200';
    %lambda_vec_string = '100:1:6000';
    lambda_vec_string = '0.0001:0.1:100';
    lambda_vec        = eval(lambda_vec_string);
    
    Ax_b_norm_1  = zeros(size(lambda_vec));
    Lx_norm_1    = zeros(size(lambda_vec));
    Ax_b_norm_2  = zeros(size(lambda_vec));
    Lx_norm_2    = zeros(size(lambda_vec));
    Ax_b_norm_3  = zeros(size(lambda_vec));
    Lx_norm_3    = zeros(size(lambda_vec));
    Ax_b_norm_4  = zeros(size(lambda_vec));
    Lx_norm_4    = zeros(size(lambda_vec));
    Ax_b_norm_5  = zeros(size(lambda_vec));
    Lx_norm_5    = zeros(size(lambda_vec));
    Ax_b_norm_6  = zeros(size(lambda_vec));
    Lx_norm_6    = zeros(size(lambda_vec));
    Ax_b_norm_7  = zeros(size(lambda_vec));
    Lx_norm_7    = zeros(size(lambda_vec));
    
    index_lambda = 1;
    
    for lambda = lambda_vec
        
        [ ridge_regression_result ]  = Regularized_Sol( Conv_Matrix        , Sim_Ct_T_noise, x0, lambda, normalize, 1);
        [ knots_coeff ]              = Regularized_Sol( Conv_Matrix*B_mat  , Sim_Ct_T_noise, knots0, lambda, normalize, 0);
        [ knots_coeff_2 ]            = Regularized_Sol( A_new_1 , Sim_Ct_T_noise, knots0, lambda, normalize, 0);
        [ knots_coeff_3 ]            = Regularized_Sol( A_new_2 , Sim_Ct_T_noise, knots0, lambda, normalize, 0);
        [ knots_coeff_PCA ]          = Regularized_Sol( A_PCA   , Sim_Ct_T_noise, knots0, lambda, normalize, 0 );
        [ knots_coeff_PCA_1 ]        = Regularized_Sol( A_PCA_1 , Sim_Ct_T_noise, knots0, lambda, normalize, 0);
        [ knots_coeff_PCA_2 ]        = Regularized_Sol( A_PCA_2 , Sim_Ct_T_noise, knots0, lambda, normalize, 0);
        
        
        % Calculate result
        b_spline_result           = B_mat        * knots_coeff;
        b_spline_result_1st_deriv = inv_matrix_1st_deriv * knots_coeff_2;
        b_spline_result_2nd_deriv = inv_matrix_2nd_deriv * knots_coeff_3;
        b_PCA_result              = B_PCA        * knots_coeff_PCA;
        b_PCA_result_1st_deriv    = inv_matrix_PCA_1st_deriv * knots_coeff_PCA_1;
        b_PCA_result_2nd_deriv    = inv_matrix_PCA_2nd_deriv * knots_coeff_PCA_2;
        
        % Remove zeros
        b_spline_result           = max(b_spline_result,0);
        b_spline_result_1st_deriv = max(b_spline_result_1st_deriv,0);
        b_spline_result_2nd_deriv = max(b_spline_result_2nd_deriv,0);
        b_PCA_result              = max(b_PCA_result,0);
        b_PCA_result_1st_deriv    = max(b_PCA_result_1st_deriv,0);
        b_PCA_result_2nd_deriv    = max(b_PCA_result_2nd_deriv,0);
        
        % Calculate norm
        Ax_b_norm_1(index_lambda) = (Conv_Matrix*ridge_regression_result - Sim_Ct_T_noise)' * (Conv_Matrix*ridge_regression_result - Sim_Ct_T_noise);
        Lx_norm_1(index_lambda)   = (ridge_regression_result)' * (ridge_regression_result);
        
        Ax_b_norm_2(index_lambda) = (Conv_Matrix*b_spline_result - Sim_Ct_T_noise)' * (Conv_Matrix*b_spline_result - Sim_Ct_T_noise);
        Lx_norm_2(index_lambda)   = (knots_coeff)' * (knots_coeff);
        
        %     Ax_b_norm_3(index_lambda) = (Conv_Matrix*b_spline_result_1st_deriv - Sim_Ct_T_noise)' * (Conv_Matrix*b_spline_result_1st_deriv - Sim_Ct_T_noise);
        %     Lx_norm_3(index_lambda)   = (deriv_matrix*b_spline_result_1st_deriv)' * (deriv_matrix*b_spline_result_1st_deriv);
        %
        %     Ax_b_norm_4(index_lambda) = (Conv_Matrix*b_spline_result_2nd_deriv - Sim_Ct_T_noise)' * (Conv_Matrix*b_spline_result_2nd_deriv - Sim_Ct_T_noise);
        %     Lx_norm_4(index_lambda)   = (deriv_matrix_2*b_spline_result_2nd_deriv)' * (deriv_matrix_2*b_spline_result_2nd_deriv);
        
        Ax_b_norm_3(index_lambda) = (Conv_Matrix*b_spline_result_1st_deriv - Sim_Ct_T_noise)' * (Conv_Matrix*b_spline_result_1st_deriv - Sim_Ct_T_noise);
        Lx_norm_3(index_lambda)   = (deriv_matrix_1st*knots_coeff_2)' * (deriv_matrix_1st*knots_coeff_2);
        
        Ax_b_norm_4(index_lambda) = (Conv_Matrix*b_spline_result_2nd_deriv - Sim_Ct_T_noise)' * (Conv_Matrix*b_spline_result_2nd_deriv - Sim_Ct_T_noise);
        Lx_norm_4(index_lambda)   = (deriv_matrix_2nd*knots_coeff_3)' * (deriv_matrix_2nd*knots_coeff_3);
        
        Ax_b_norm_5(index_lambda) = (Conv_Matrix*b_PCA_result - Sim_Ct_T_noise)' * (Conv_Matrix*b_PCA_result - Sim_Ct_T_noise);
        Lx_norm_5(index_lambda)   = (knots_coeff_PCA)' * (knots_coeff_PCA);
        
        Ax_b_norm_6(index_lambda) = (Conv_Matrix*b_PCA_result_1st_deriv - Sim_Ct_T_noise)' * (Conv_Matrix*b_PCA_result_1st_deriv - Sim_Ct_T_noise);
        Lx_norm_6(index_lambda)   = (PCA_deriv_matrix_1*knots_coeff_PCA_2)' * (PCA_deriv_matrix_1*knots_coeff_PCA_2);
        
        Ax_b_norm_7(index_lambda) = (Conv_Matrix*b_PCA_result_2nd_deriv - Sim_Ct_T_noise)' * (Conv_Matrix*b_PCA_result_2nd_deriv - Sim_Ct_T_noise);
        Lx_norm_7(index_lambda)   = (PCA_deriv_matrix_2*knots_coeff_PCA_2)' * (PCA_deriv_matrix_2*knots_coeff_PCA_2);
        
        index_lambda = index_lambda + 1;
        
    end
    
    % Plot L curve
    fig_num = figure;
    subplot(7,1,1);
    h1 = loglog(Ax_b_norm_1,Lx_norm_1);
    xlabel('||AR-C||2');
    ylabel('||R||2');
    title(['L-Curve - Simple regularization. Lambda:' lambda_vec_string],'FontWeight','bold');
    
    subplot(7,1,2);
    h1 = loglog(Ax_b_norm_2,Lx_norm_2);
    xlabel('||ABV-C||2');
    ylabel('||V||2');
    title(['L-Curve - Spline regularizatio. Lambda:' lambda_vec_string],'FontWeight','bold');
    
    subplot(7,1,3);
    h1 = loglog(Ax_b_norm_3,Lx_norm_3);
    xlabel('||ABV-C||2');
    ylabel('||LV||2');
    title(['L-Curve - Spline 1st deriv. regularization. Lambda:' lambda_vec_string],'FontWeight','bold');
    
    subplot(7,1,4);
    h1 = loglog(Ax_b_norm_4,Lx_norm_4);
    xlabel('||ABV-C||2');
    ylabel('||LLV||2');
    title(['L-Curve - Spline 2nd deriv. regularization. Lambda:' lambda_vec_string],'FontWeight','bold');
    
    subplot(7,1,5);
    h1 = loglog(Ax_b_norm_5,Lx_norm_5);
    xlabel('||ABV-C||2');
    ylabel('||LLV||2');
    title(['L-Curve - PCA regularization. Lambda:' lambda_vec_string],'FontWeight','bold');
    
    subplot(7,1,6);
    h1 = loglog(Ax_b_norm_6,Lx_norm_6);
    xlabel('||ABV-C||2');
    ylabel('||LLV||2');
    title(['L-Curve - PCA 1st deriv. regularization. Lambda:' lambda_vec_string],'FontWeight','bold');
    
    subplot(7,1,7);
    h1 = loglog(Ax_b_norm_7,Lx_norm_7);
    xlabel('||ABV-C||2');
    ylabel('||LLV||2');
    title(['L-Curve - PCA 2nd deriv. regularization. Lambda:' lambda_vec_string],'FontWeight','bold');
    
    
    % Print result to PDF
    Fig_name = [filter_type '_L_Curves.png'];
    Title    = [filter_type ' L-Curves'];
    subTitle = '';
    [idx_fig] = Print2Pdf(fig_num, idx_fig, Fig_name, './Run_Output/',Title, subTitle);
    
    % Find curvature point (should do this manually)
    
    % L1
    curvature_y_val = 184.6;
    [value idx ] = min(abs(Lx_norm_1 - curvature_y_val));
    needed_lambda = lambda_vec(idx);
    % L2
    curvature_y_val = 120.4;
    [value idx ] = min(abs(Lx_norm_2 - curvature_y_val));
    needed_lambda = lambda_vec(idx);
    % L3
    curvature_y_val = 117.5;
    [value idx ] = min(abs(Lx_norm_3 - curvature_y_val));
    needed_lambda = lambda_vec(idx);
    % L4
    curvature_y_val = 0.2588;
    [value idx ]    = min(abs(Lx_norm_4 - curvature_y_val));
    needed_lambda   = lambda_vec(idx);
    % L5
    curvature_y_val = 170.1;
    [value idx ]    = min(abs(Lx_norm_5 - curvature_y_val));
    needed_lambda   = lambda_vec(idx);
    % L6
    curvature_y_val = 26.86; % Or 585.4
    [value idx ]    = min(abs(Lx_norm_6 - curvature_y_val));
    needed_lambda   = lambda_vec(idx);
    % L7
    curvature_y_val = 217.6;
    [value idx ]    = min(abs(Lx_norm_7 - curvature_y_val));
    needed_lambda   = lambda_vec(idx);
    
end

% Assign return values
if RealData_Flag
    ridge_regression_result    = NaN;
    b_spline_result            = NaN;
    b_spline_result_1st_deriv  = NaN; 
    %b_spline_result_2nd_deriv  = NaN; 
    b_PCA_larss_result         = NaN;
    b_PCA_result_1st_deriv     = NaN; 
    b_PCA_result_2nd_deriv     = NaN; 
    %idx_fig                    = NaN;
end



end

