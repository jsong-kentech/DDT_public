% This code calculates a distribution of diffusion times (DDT) for a given
% EIS data. Written by Juhyun Song, 2018.
% Reference: Juhyun Song and Martin Z. Bazant, Phys. Rev. Lett. 120, 116001

% noCV: not include cross-validation to find optimal lambda and R_intercept

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
        R_star = 4.46; % real-axis intercept of diffusion impedance, predetermined [ohm]
        l_star = -1;  % regularization parameter (lambda), log scale, predetermined    
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





%% DDT Inversion

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






 %% DDT Inversion: Confidence interval - Bootstrap method

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



