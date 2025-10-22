function [K,D] = DDT_public_Kernel(type_kernel,x,ppd_t,t)

    M =length(x); % length of the experimental data vectors
    dt=1/ppd_t;
    N=length(t);
    
    % K_matrix
    K = zeros(M,N);
    if type_kernel ==1  % Finite Warburg kernel (transmissive boundary), Planar
        for m =1:M
            for n=1:N
                K(m,n) = dt*sqrt(1i*exp(t(n)-x(m)))*(tanh(sqrt(1i*exp(t(n)-x(m)))));
            end
        end
    elseif type_kernel == 2 % Bounded Warburg kernel (blocking boundary), Cylindrical
        for m =1:M
            for n=1:N
                if (t(n)-x(m)) > 6
                % K_D(m,n) = dtD*sqrt(1i*exp(tD(n)-x(m)))*(tanh(sqrt(1i*exp(tD(n)-x(m))))); 
                K(m,n) = dt*((1/sqrt(2*exp(t(n)-x(m))))*(1-1i)-(1/(2*exp(t(n)-x(m))))*1i)^-1;                
                else
                K(m,n) = dt*sqrt(1i*exp(t(n)-x(m)))*(besseli(1,(sqrt(1i*exp(t(n)-x(m)))))/besseli(0,(sqrt(1i*exp(t(n)-x(m))))));     
                end
            end
        end
    elseif type_kernel == 3 % Finite Warburg kernel (transmissive boundary), Planar
        for m =1:M
            for n=1:N
                K(m,n) = dt*sqrt(1i*exp(t(n)-x(m)))*coth(sqrt(1i*exp(t(n)-x(m)))); 
            end
        end
    end
    

    
    % D_matrix
    
    D = zeros(N,N);
    for n = 1:N
        if n ==1 % forward
        D(n,n)=1;
        D(n,n+1)= -2;
        D(n,n+2)=1;
        elseif n == N % backward
        D(n,n-2)=1;
        D(n,n-1)= -2;
        D(n,n)=1;    
        else % central
        D(n,n-1)=1;
        D(n,n)= -2;
        D(n,n+1)=1;
        end
    end
    
    D = D/(dt);



end