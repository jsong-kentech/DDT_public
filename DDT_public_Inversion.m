function [r] = DDT_public_Inversion(y,t,K,W_re,W_im,D,...
                            l,...
                            scaling_real,scaling_imag)
    
    lambda=10^(l);
    y_RE = real(y); y_IM = imag(y);
    N=length(t);
    K_RE = real(K);
    K_IM = imag(K);

    % DEFINE THE QUAD PROGRAMMING

    H = scaling_real*K_RE'*(W_re'*W_re)*K_RE + scaling_imag*K_IM'*(W_im'*W_im)*K_IM + lambda*(D'*D);
    f = - scaling_real*y_RE'*(W_re'*W_re)*K_RE - scaling_imag*y_IM'*(W_im'*W_im)*K_IM; 

    options = optimoptions('quadprog','algorithm','interior-point-convex',...
    'maxiterations',1000,'optimalitytolerance',eps,'steptolerance',eps,...
    'constrainttolerance',eps);

    lb=zeros(N,1);
    % q = quadprog(H,f,[],[],[],[],[],[],[],options); % without positive constraint
    r = quadprog(H,f,[],[],[],[],lb,[],[],options); % with positive constraint.


end