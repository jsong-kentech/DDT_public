function [W_re,W_im] = DDT_public_Weight(scale,M,type_weight)
    
    W_re = zeros(M);
    W_im = zeros(M);
    if type_weight == 0
        for m = 1:M
            W_re(m,m) = 1; % only diagonal components
        end
        W_im = W_re;
    elseif type_weight ==1
        for m = 1:M
        W_re(m,m) = abs((scale(m)).^-1); % only diagonal components
        end
        W_im = W_re;
    elseif type_weight ==2
        for m = 1:M
        W_re(m,m) = abs((real(scale(m))).^-1); % only diagonal components
        W_im(m,m) = abs((imag(scale(m))).^-1); 
        end  
    end
end

