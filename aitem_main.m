function [mom_values,b_values,wave_speeds,dcdp,...
    growth_rates,err,solutions] = aitem_main(n,xl,n_mom_values,n_b_values_side,pmax,bmax,c,...
    lambda_1,lambda_2,sigma1,sigma2,mu1,mu2,nu11,nu12,nu21,nu22,gamma1,gamma2,hmin)


errcalc = 1;         		% calculate error of solution

%
%   renormalisation
%

alfa = 1./sqrt(abs(caitem_in.mu1*caitem_in.mu2));
beta = 1./sqrt(abs(caitem_in.mu1*caitem_in.mu2));

sigma1 = sigma1*beta/alfa;
sigma2 = sigma2*alfa/beta;
mu1 = mu1*alfa;
mu2 = mu2*beta;
nu11 = nu11/alfa*beta^2;
nu21 = nu21*beta;
nu12 = nu12*alfa;
nu22 = nu22/beta*alfa^2;
gamma1 = gamma1/alfa;
gamma2 = gamma2/beta;

%
% do a check on eigs here, if negative: swap coefficients
% if one positive and one negative: bail
%
% find eigenvalues of M
%

M = [lambda_1 sigma1; 
     sigma2   lambda_2];

detM = det(M);

assert(detM>0,'determinant of M is not positive!');
assert(isreal(M), 'eigenvalues of M are not real!');

if (trace(M) < .0)
    c = -c; sigma1 = -sigma1; sigma2 = -sigma2; lambda_1 = -lambda_1; lambda_2 = -lambda_2;
    mu1 = -mu1; nu11 = -nu11; nu21 = -nu21; gamma1 = -gamma1;
    mu2 = -mu2; nu22 = -nu22; nu12 = -nu12; gamma2 = -gamma2;
end

M = [lambda_1 sigma1; sigma2 lambda_2];
NL = @(A,B,b) [ c*A+mu1*A.^2+nu11*B.^2+2*nu21*A.*B+gamma1*b; ...
               -c*B+mu2*B.^2+nu22*A.^2+2*nu12*A.*B+gamma2*b];

dx = xl/n;
x = dx*(1:n);
uo = (1-(tanh(x-.5*xl)).^2);  % initial solution is this times b


% 
% initialize the Fourier differentiation matrix
%

ind = 1:n-1;
r = [0 .5*(-1).^(ind+1)./tan(pi*ind/n)];
col = [0 r(n:-1:2)];
Df = 2*pi/xl*toeplitz(col,r);
kk = complex(0,1)*2*pi/xl*[0:(n/2-1) (-n/2):(-1)];
k2 = kk.^2;

%
% work from bo = 0 out
% need some estimate of appropriate size for bmax and pm
%
% to look for unforced solitary waves on b0=0, set n_b_values_side = 0 and bmax = 0
%

% suppose pmax=1, n_mom_values=100; then have grid [0.01, 1] with 100
% points (0 is NOT included):
mom_values = linspace(pmax/n_mom_values,pmax,n_mom_values);

n_b_values = 2*n_b_values_side+1;

% b always in symmetric interval [-bmax, bmax]
b_values = linspace(-bmax,bmax,n_b_values);

% there will be 1 wavespeed for every momentum and b
wave_speeds = nan(n_mom_values,n_b_values);

hs = nan(n_mom_values,n_b_values);

% will store solutions for every momentum and b in a multidimensional array
solutions = nan(n_mom_values,n_b_values,2,n);

A0 = unifrnd(-1,1);
B0 = unifrnd(-1,1);
U0 = [A0*uo; B0*uo];

h = hmin; % start with time step hmin

% start going through all momentum values
for idx_p = 1:n_mom_values
    
    P = mom_values(idx_p);
    
    idx_b = n_b_values_side+1; % this is index where b=0
    
    b = b_values(idx_b)*uo; % this is what topography is
    
    [hf,nit,wave_speeds(idx_p,idx_b),solution_by_aitem] = caitem(P,b,dx,h,M,NL,U0);
	
    if (wave_speeds(idx_p,idx_b) > -c)
		wave_speeds(idx_p,idx_b) = nan;
    end
     
    % rerun with different initial amplitudes
    
    kk = .0;  % number of repeats
    
    while (isnan(wave_speeds(idx_p,idx_b)) && kk < 4)
        
        h = hmin;
	    A0 = unifrnd(-1,1);
        B0 = unifrnd(-1,1); 
        
        U0 = [A0*uo; B0*uo];
        
        [hf,nit,wave_speeds(idx_p,idx_b),solution_by_aitem] = caitem(P,b,dx,h,M,NL,U0);
        
		if (wave_speeds(idx_p,idx_b) > -c)
			wave_speeds(idx_p,idx_b) = nan;
        end
        
        kk = kk+1;
        
    end
    
    % done with rerunning;
    % if the wavespeeds are still either >-1 or NaNs, make the solution
    % NaNs (i.e. we admit we couldn't find it)
    
    if (isnan(wave_speeds(idx_p,idx_b)))
        solutions(idx_p,idx_b,:,:) = nan(size(solution_by_aitem));
    end
    
    % save time step for the current values of momentum and b
    hs(idx_p,idx_b) = hf;
    
    h = hf;
    
    U0 = solution_by_aitem;
    
    % save the solution found by AITEM
    solutions(idx_p,idx_b,:,:) = solution_by_aitem;

end  % stop going through momentum values; recall this was only for b=0!

% start going through the values of momentum again

for idx_p = 1:n_mom_values
    
    P = mom_values(idx_p);
    
    h = hs(idx_p,idx_b); % idx_b initially corresp. to b=0
    
    U0 = squeeze(solutions(idx_p,idx_b,:,:)); % start from the already found solution
%
% then work out on each side
%
    for ia = 1:2
        
        hu = h; % time step to supply to aitem
        solution_to_start_aitem = U0;
        
        for jb = 1:n_b_values_side
            
            % while ia=1 index ib is n_b_values_side+1-jb, which is 
            % n_b_values_side, n_b_values_side-1,.., 1 - going from b=0
            % down to the very first index; once a1=2, start going up to
            % 2*n_b_values_side+1 (the very top index)
            
            ib = n_b_values_side+1+(-1)^(ia)*jb; 
            
            b = b_values(ib)*uo;
            
            [hf,nit,wave_speeds(idx_p,ib),solution_by_aitem] = caitem(P,b,dx,hu,M,NL,solution_to_start_aitem);
			
            if (wave_speeds(idx_p,ib) > -c)
				wave_speeds(idx_p,ib) = nan;
            end
            
            % like before, try repetitions
            
            kk = .0;
            
            while (isnan(wave_speeds(idx_p,ib)) && kk < 4)

                h = hmin;

				A0 = unifrnd(-1,1);
                B0 = unifrnd(-1,1);

                U0 = [A0*uo; B0*uo];
                
                [hf,nit,wave_speeds(idx_p,ib),solution_by_aitem] = caitem(P,b,dx,h,M,NL,U0);
                
                if wave_speeds(idx_p,ib)>-c
                	wave_speeds(idx_p,ib)=nan;
                end
                
				kk = kk+1;
            end
            
            
            if (isnan(wave_speeds(idx_p,ib)))
                solutions(idx_p,ib,:,:) = nan(size(solution_by_aitem));
                break; % get out of looping on this side and got to the other size (change ai)
            end
            
            % however, if solution has been found, save time step
            hs(idx_p,ib) = hf;
            hu = hf;
            solution_to_start_aitem = solution_by_aitem; % solution to supply to aitem
            
            solutions(idx_p,ib,:,:) = solution_by_aitem;
        end
    end
end

%
% try to fill in any NaN calculations
%

for ib = 1:n_b_values
    
    b = b_values(ib)*uo;
    
    for idx_p = 2:n_mom_values
        
        P = mom_values(idx_p);
        
        if (isnan(wave_speeds(idx_p,ib)) && ~isnan(wave_speeds(idx_p-1,ib)))
            
            U0 = squeeze(solutions(idx_p-1,ib,:,:));
            
            h = hs(idx_p-1,ib);
            
            [hs(idx_p,ib),nit,wave_speeds(idx_p,ib),solution_by_aitem] = caitem(P,b,dx,h,M,NL,U0);
            
            solutions(idx_p,ib,:,:) = solution_by_aitem;
        end
    end
end

%
% calculate derivative of c wrt P (should give instability criteria)
%

dcdp = zeros(size(wave_speeds));

for i = 1:length(b_values)
    dcdp(2:n_mom_values-1,i) = (wave_speeds(3:n_mom_values,i)-wave_speeds(1:n_mom_values-2,i))./(mom_values(3:n_mom_values)-mom_values(1:n_mom_values-2))'; 
end

%
% calculate the norm of the error and the growthrate
%

growth_rates = zeros(n_mom_values,n_b_values);

err = zeros(n_mom_values,n_b_values);

for idx_p = 1:n_mom_values
    
    for ib = 1:n_b_values
        
        if isnan(wave_speeds(idx_p,ib))
            growth_rates(idx_p,ib) = NaN;
            err(idx_p,ib) = NaN;
        else
            as = squeeze(solutions(idx_p,ib,1,:));
            bs = squeeze(solutions(idx_p,ib,2,:));
            Delta = wave_speeds(idx_p,ib);
            
            %
            %   calculate norm of error
            %
            if (errcalc)
                
                asxx = real(ifft(k2'.*fft(as))); % Axx
                bsxx = real(ifft(k2'.*fft(bs))); % Bxx
                
                E1 = (Delta+c)*as+(lambda_1*asxx+sigma1*bsxx)+mu1*as.^2+nu11*bs.^2 ...
                    +2*nu21*as.*bs+gamma1*b_values(ib)*uo';
                
                E2 = (Delta-c)*bs+(lambda_2*bsxx+sigma2*asxx)+mu2*bs.^2+nu22*as.^2 ...
                    +2*nu12*as.*bs+gamma2*b_values(ib)*uo';
                
                err(idx_p,ib) = norm([E1' E2']);  % returns maximum singular value (i.e. 2-norm)
                
            end
            %
            %   calculate growthrate
            %						
				
            M11 = Df*(diag(Delta+c+2*mu1*as+2*nu21*bs)+lambda_1*Df^2);
            M12 = Df*(diag(2*nu11*bs+2*nu21*as)+sigma1*Df^2);
            M21 = Df*(diag(2*nu12*bs+2*nu22*as)+sigma2*Df^2);
            M22 = Df*(diag(Delta-c+2*mu2*bs+2*nu12*as)+lambda_2*Df^2);
            
            Mt = -[M11 M12; M21 M22];
            
            meigs = eig(Mt);

            % store not just the maximal real part but the whole
            % eigenvalue with the maximal real part

            [sig, idx] = max(real(meigs));
            growth_rates(idx_p,ib) = meigs(idx);


        end  % if isnan(wave_speeds(idx_p,ib))...
    end    
end

    
function [dt,ic,wave_speed1,U] = caitem(P,xf,dx,dtc,m,g,Uo)
%
% function to do equivalent to AITEM method for coupled KdV equations:
%
% idea is to do implicit (A_xx)-explicit (-wave_speed1*A+g(A)) differencing and 
% choose optimal value of delta t
%
    hmax = 10.;
    sf = 1.05;
    nt = 1000;
    conv = 1.e-5;
    rmsp = 10.;
%
    N = length(xf);
%
    T1 = zeros(size(xf));
     
    if N > 1 
        L = N*dx;
    else
        disp('Size of x field needs to be larger than one');
        return
    end
%
% initialize wavenumbers and inversion matrices
%
    km = 2*pi/L;
    K = km*[0:(N/2-1) (-N/2):(-1)];
    Lo = -K.^2;
%
% initialize field if Uo is not empty
%
    if isempty(Uo)
        disp('Uo is empty');
        return;
    else
        for j = 1:N
            T1(j) = Uo(1,j)^2+Uo(2,j)^2;
        end
        Po = dx*sum(T1);
        U = sqrt(P/Po)*Uo;
    end
%
% start of iteration routine, continue until |dU| < convcrit
%
    dt = max(dtc,hmin);
    
    for ic = 1:nt
        
        Up = U;
%
%   use dt
%
        [wave_speed1,dPdt,rmse,U] = caitem_update(dt,dx,N,P,xf,Up,g,m,Lo);
        if (rmse <= conv)
          break
        elseif (rmse >= rmsp)
            dt = max(dt/sf^2,hmin);
        else
%
%   use sf*dt
%
        [wave_speed2,dPdt,RMSE_ERROR,U2] = caitem_update(sf*dt,dx,N,P,xf,Up,g,m,Lo);
%
% following criteria seems to be faster than dP2 > dP1, and there is
% evidence to support this
%
            if (wave_speed2 < wave_speed1)
                
                U = U2;
                
                wave_speed1 = wave_speed2; % make it equal to smaller wave speed
                
                dt = min(sf*dt,hmax);
                
                rmse = RMSE_ERROR;
            end
        end
        
        wave_speed1 = NaN;
        rmsp = rmse;
        
        if (dt >= hmax)
            return;
        end
    end
end

end
