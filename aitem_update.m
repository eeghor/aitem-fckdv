%%%% 
%%%% AITEM UPDATE for a system of forced cKdV equations
%%%%

function [wave_speed,dPdt,RMSE_ERROR,Un] = aitem_update(c,mu1,mu2,nu11,nu21,nu22,nu12,...
                            gamma1,gamma2,dt,dx,N,Momentum_Uo,bo,Uo,M,mK2)
 
        % M = [lambda1 sigma1; sigma2 lambda2];
        % Uo is an array like row 1: [A], row 2: [B]
        % c is normally 1
        % mK2 is -k.^2
        
        A = Uo(1,:); B = Uo(2,:);
        
        L00 = zeros(size(Uo)); 
        T1 = zeros(size(Uo)); 
        T2 = zeros(size(Uo));
        
        Nonl_Part = [ c*A+mu1*A.^2+nu11*B.^2+2*nu21*A.*B+gamma1*bo; ...
                     -c*B+mu2*B.^2+nu22*A.^2+2*nu12*A.*B+gamma2*bo];
        
        FFT_Nonl_Part = fft(Nonl_Part,[],2);
        
        FFT_Uo = fft(Uo,[],2);
        
        for j = 1:N
            L00(:,j) = mK2(j)*M*FFT_Uo(:,j)+FFT_Nonl_Part(:,j);
        end
%
% calculate M^{-1}*U and M^{-1}*U_t
%
        for j = 1:N
            Minv = (eye(2)-dt*mK2(j)*M)\eye(2)*dt;
            T1(:,j) = Minv*L00(:,j);
            T2(:,j) = Minv*FFT_Uo(:,j);
        end
        
        MUt = real(ifft(T1,[],2));
        MU = real(ifft(T2,[],2));
%
% calculate wavespeed
%
        wave_speed = -sum(Uo(1,:).*MUt(1,:)+Uo(2,:).*MUt(2,:)) ...
              /sum(Uo(1,:).*MU(1,:)+Uo(2,:).*MU(2,:));
%
% calculate update
%
        Us = Uo+MUt+wave_speed*MU;
%
% normalize
%
        Momentum_Us = dx*sum(Us(1,:).^2+Us(2,:).^2);
        Un = sqrt(Momentum_Uo/Momentum_Us)*Us;
        dPdt = (Momentum_Us-Momentum_Uo)/dt; 
%
% calculate rms error
%
        dut = (Un-Uo)/dt;
        
        for j = 1:N
           A(j) = dut(1,j)^2+dut(2,j)^2;
        end
        
        RMSE_ERROR = sqrt(dx*sum(A));
end