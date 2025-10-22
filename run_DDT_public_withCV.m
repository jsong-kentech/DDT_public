% This code calculates a distribution of diffusion times (DDT) for a given
% EIS data. Written by Juhyun Song, 2018.
% Reference: Juhyun Song and Martin Z. Bazant, Phys. Rev. Lett. 120, 116001

% withCV: include cross-validation to find optimal lambda and R_intercept


clear; clc; close all
%% Configuration

% Data

    % Data path and name
    filename_data = 'SiNW_data.mat';
    
    % Load data (see "code_to_load_data_and_hyperparameters_other_concentrations.txt" for other data set)
    data = load(filename_data);    
        w_data= data.w_1274; % frequency data (vector) [Rad/sec]
        z_data = data.z_1274; % impedance data (vector, complex) [Ohm]
        n_begin = 40; % index of the first data point to include in DDT analysis (only include the diffusion part)
        n_end = 57; % index of the last data point to include 

   % hyperparameters
        R_lb = 4.4; % real-axis intercept of diffusion impedance, lower bound in cross-validation [ohm]
        R_ub = 4.5; % real-axis intercept of diffusion impedance, uppder bound in cross-validation [ohm]
        l_lb = -2;  % regularization parameter (lambda), log scale, lower bound in cross-validation 
        l_ub = -0.5;% regularization parameter (lambda), log scale, upper bound in cross-validation    
        t_lb = 0; % time scale to include in DDT analysis, lower limit, log scale
        t_ub = 5; % time scale to include in DDT analysis, upper limit, log scale  


%  Config DDT

    % Type of kernel function
    type_kernel = 2;    % 1, BD (Bounded diffusion) Planar 
                        % 2, BD (Bounded Diffusion) Cylindrical
                        % 3, FLW (transmissive diffusion) Planar
    % Type of weighting 
    type_weight = 1; % 0, Uniform / 1, Relative 

    
%% Data Pre-processing
    
    % plot the imported impedance data
    figure(1);
    hold on;
    plot(real(z_data),-imag(z_data),'ko')
    axis([0,7,0,7])
    xlabel('Z_{re} [Ohm]')
    ylabel('-Z_{im} [Ohm]')
    legend_fig1{1} = 'exp data';
    legend(legend_fig1,'location','northwest')

    % Admittance 
    y_data = z_data.^-1;
    
    % Separate diffusion part
    z_D = z_data(n_begin:n_end,1); 
    y_D = z_D.^-1;
    w_D = w_data(n_begin:n_end,1);
        % relaxation data are not used


%% Preparing the Inversion

    % Define the tD and tauD 
        % point per decade in tD vector
        ppd_t = (n_end - n_begin + 5)/(t_ub - t_lb); % added 6 more points to make the inversion overdetermined
        % generate tD vector and tauD vector
        t = t_lb:1/ppd_t:t_ub;
        tau = exp(t); % tau is define in Reference page 2, t = log(tau)
        N = length(t);

    % Calculate the Kernels
        % frequency and x
        x = -log(w_D); % x is defined in Reference page 2, x = -log(w)
        M = length(x);

        % Kernel
        [K,D] = DDT_public_Kernel(type_kernel,x,ppd_t,t);



%% DDT Inversion: (1) Course Optimization of Hyperparameters (Cross Validation)

    % Hyperparameter grid to evaluate
        % lambda (regularization parameter) grid
        l_grid = linspace(l_lb,l_ub,20);
        L_l = length(l_grid);
        % R (intercept resistance) grid
        R_grid = linspace(R_lb,R_ub,20);
        L_R = length(R_grid);

    % Evaluate cross-validation error on the grid
        CVE_grid = zeros(L_l, L_R);
        for i_l = 1:L_l
            for i_R = 1:L_R

                CVE_grid(i_l,i_R) = DDT_public_CrossValidation(y_D,t,K,D,...
                     [l_grid(i_l);R_grid(i_R)],M,type_weight);        


            end
        end

    % Look up the course optimum
       [CVE_hat,ind_hat] = min(CVE_grid(:));
       [row_hat,col_hat] = ind2sub(size(CVE_grid),ind_hat);
       l_hat = l_grid(row_hat);
       R_hat = R_grid(col_hat);

    % Plot the course minimum
        figure(2); hold on;
        surf(R_grid,l_grid,CVE_grid)
        plot3(R_hat,l_hat,CVE_hat,'ok','MarkerFaceColor','k')
        view(-35,30)
        axis([R_lb,R_ub,l_lb,l_ub,0.9*min(min(CVE_grid)),1.1*max(max(CVE_grid))])
        xlabel('R [Ohm]')
        ylabel('l=log(\lambda)')
        zlabel('CVE')
        legend_fig2 ={'CVE_{grid}','CVE_{hat}'};
        legend(legend_fig2)
       
    
        figure(3); hold on;
        plot(l_grid,CVE_grid(:,col_hat),'-b')
        plot(l_hat,CVE_hat,'ob','MarkerFaceColor','b')
        axis([l_lb,l_ub,0.9*min(CVE_grid(:,col_hat)),1.1*max(CVE_grid(:,col_hat))])
        xlabel('l=log(\lambda)')
        ylabel('CVE')
        legend_fig3 = {'CVE(R=R_{hat})','CVE_{hat}'};
        legend(legend_fig3)


 %% DDT Inversion: (2) Fine Optimization of Hyperparameters (Cross Validation)

    % Define the cost function and the variables to optimize by function
    % handle
        CVE_fhandle = @(para)DDT_public_CrossValidation(y_D,t,K,D,para,M,type_weight);
                        % para is 2-comp vector [l,R]
    % Optimization options
        CVE_option = optimoptions('fminunc','OptimalityTolerance',eps,'StepTolerance',eps,'FunctionTolerance',eps);
    % Initial guess
        l_0 = l_hat;
        R_0 = R_hat;
    % Run optimization
        [para_star,CVE_star] = fminunc(CVE_fhandle,[l_0,R_0],CVE_option);
        l_star = para_star(1);
        R_star = para_star(2);

    % Plot the fine minimum
        figure(2); hold on;
        plot3(R_star,l_star,CVE_star,'or','MarkerFaceColor','r')
        legend_fig2{end+1} = 'CVE_{star}';
        legend(legend_fig2)




%% DDT Inversion: (3) Inversion with the optimal hyperparameters

    % Pre-processing: horizontal shift of diffusion EIS data
        y_D_adj = (z_D - R_star).^-1;
    % Calculate weighting matrix
        [W_re,W_im] = DDT_public_Weight(y_D_adj,M,type_weight);
    % Inversion by quadratic programming
        q_star = DDT_public_Inversion(y_D_adj,t,K,W_re,W_im,D,l_star,1,1);
    % Fitted model evaluation
        y_D_adj_star = K*q_star;
        z_D_star = y_D_adj_star.^-1 + R_star;
    % Normalization
        q_star = q_star./sum(q_star)*ppd_t;


    % Plot the solution
        figure(1); hold on;
        plot(real(z_D_star),-imag(z_D_star),'ob');
        legend_fig1{end+1} = 'DDT model';
        legend(legend_fig1)

        figure(4); hold on;
        plot(t,q_star,'ob')
        axis([t_lb,t_ub,0,1.1*max(q_star)]);
        xlabel ('t = log(\tau)')
        ylabel ('q = \tauP(\tau)')
        legend_fig4{1} = 'DDT solution';
        legend(legend_fig4)






 %% DDT Inversion: (4) Confidence interval - Bootstrap method

    % prepare Bootstrap method
        N_resample = 200; % number of re-sample (multiple of 5)
        M_resample = M; % number of data points in each sample = number of original data
        idx_resample = unidrnd(M_resample,[M_resample,N_resample]);
                % N columns, each column index of resampling (M-rows), random numbers from 1 - M

    % Resampling
        q_resample = zeros(N,N_resample);
        for i_resample = 1:N_resample
            y_resample = y_D_adj(idx_resample(:,i_resample)); % resampled data
            W_re_resample = W_re(idx_resample(:,i_resample),idx_resample(:,i_resample));
            W_im_resample = W_im(idx_resample(:,i_resample),idx_resample(:,i_resample));
            K_resample = K(idx_resample(:,i_resample),:); % kernel matrix for resample

            % DDT solution of resampled data
            q_resample(:,i_resample) = DDT_public_Inversion(y_resample,t,K_resample,...
                W_re_resample,W_im_resample,D,l_star,1,1);

            % Normalization
            q_resample(:,i_resample) = q_resample(:,i_resample)./sum(q_resample(:,i_resample))*ppd_t;

        end

    % Plot resampled solutions, superposition 
    figure(5); hold on;
    for i_resample = 1:N_resample
        plot(t,q_resample(:,i_resample))
    end
    axis([t_lb,t_ub,0,1.1*max(max(q_resample))]);
    xlabel ('t = log(\tau)')
    ylabel ('q = \tauP(\tau)')

    % Rank resampled solutions
    q_resample_sort = zeros(size(q_resample));
    q_top5p = zeros(N,1);
    q_bot5p = zeros(N,1);
     for n = 1:N % for each t, sorting 
         [q_resample_sort(n,:)] = sort(q_resample(n,:)); % sorting all resampling at spedicif t
         % top 5%
            k_top5p = round(N_resample*0.95);
            q_top5p(n) = q_resample_sort(n,k_top5p);
         % bottom 5%
            k_bot5p = round(N_resample*0.05);
            q_bot5p(n) = q_resample_sort(n,k_bot5p);
     end

    % Error bar
    q_ebar_pos = q_top5p - q_star;
    q_ebar_neg = q_star - q_bot5p;

    figure(4)
    errorbar(t,q_star,q_ebar_neg,q_ebar_pos,'sq')
    axis([t_lb,t_ub,0,1.1*max(q_top5p)]);
    legend_fig4{end+1} = 'Error bar';
    legend(legend_fig4)


% end.



