function [CVE,r_real,r_imag] = DDT_public_CrossValidation(y_exp,t,K,D,...
        para,M,type_weight)

    
    % v2
    l_now = para(1);
    R_now = para(2);
    z_adj = y_exp.^-1 - R_now;
    y_exp = z_adj.^-1; % data shifted by the intercept resistance

    % weighting matx
    scale_vector = y_exp;
    [W_re,W_im] = DDT_public_Weight(scale_vector,M,type_weight);
    
    %% INVERSION ( y -> Q ) BASED ON REAL
    % prescribe the t-values on which the inverse distribution is evaluated
    [r_real] = DDT_public_Inversion(y_exp,t,K,W_re,W_im,D,...
                                l_now,1,0);
    % calculate CVE_imag, using q_real
    CVE_imag = norm(W_im*(imag(y_exp) - imag(K*r_real)));
    
    %% INVERSION BASED ON IMAGINARY
    [r_imag] = DDT_public_Inversion(y_exp,t,K,W_re,W_im,D,...
                                l_now,0,1);
    % calculate CVE_real, using q_imag
    CVE_real = norm(W_re*(real(y_exp) - real(K*r_imag)));
    
    %% Total CVE
    CVE  = CVE_imag + CVE_real;
        
end