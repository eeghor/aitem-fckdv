%%%% 
%%%% FIND KISSING POINTS and EVALUATE SYSTEM COEFFICIENTS at these points
%%%%

close all % close all figures (if any)

%%% PART I: FIND A KISSING POINT

% colors to be used for plots
blewish=[13 140 204]./255;
darkgreen=[23 135 0]./255;
greybrown=[103 101 110]./255;
orange=[223 138 74]./255;
blu=[20 68 246]./255;
lightgray=[0.9 0.9 0.9];

fontsz=14;
MESH_DPS = 300;  % mesh density (points)

% put the model parameter values in a structure
s=struct('U0', 0, ...
         'U1', -0.2,...
         'U2', 0, ...
         'U3', 0, ...
         'H1', 0.1, ...
         'H2', 0.6,...
         'H3', 0.3, ...
         'H12', 0, ...
         'H32', 0,...
         'g12', 0.6,...
         'g32', 0.4);
     
s.('H12')=s.('H1')/s.('H2');
s.('H32')=s.('H3')/s.('H2');

% range of U3 values
U3from=-1;
U3to=1;
U3step=1e-4;
U3rng=U3from:U3step:U3to;

% check if these values make sense
assert(s.('H1')>0, 'H1 MUST be positive!')
assert(s.('H2')>0, 'H2 MUST be positive!')
assert(s.('H3')>0, 'H3 MUST be positive!')
assert(s.('g12')>=0, 'g12 MUST be non-negative!')
assert(s.('g32')>=0, 'g32 MUST be non-negative!')
assert(s.('H1')+s.('H2')+s.('H3')==1, 'H1+H2+H3 MUST be 1!')

% check if g12+g32=1; if not, do rescaling
if s.('g12')+s.('g32') ~= 1  
    gsk=sqrt(s.('g12')+s.('g32'));
    s.('g12')=s.('g12')/gsk^2;
    s.('g32')=s.('g32')/gsk^2;
    s.('U0')=s.('U0')/gsk;
    s.('U1')=s.('U1')/gsk;
    s.('U2')=s.('U2')/gsk;
    s.('U3')=s.('U3')/gsk;
    U3rng=U3rng/gsk;
    fprintf(1,'done rescaling to ensure that g12+g32=1...\n');
end

% set up the z axis: it's [-H1,H2+H3] which has the length 1
dz=1e-4;  % length of a single interval 
nc=round(1/dz+1);  % number of points incl. boundary points
z=linspace(-s.('H1'),s.('H2')+s.('H3'), nc); 
zlow=(z<0);  % boolean vector with 1s only where z in lower layer
zmid=(z>=0 & z<=s.('H2'));
ztop=(z>s.('H2'));

% create the base velocity profile
Omega1=@()   (s.('U1')-s.('U0'))/s.('H1');
Omega2=@()   (s.('U2')-s.('U1'))/s.('H2');
Omega3=@(u3) (u3-s.('U2'))/s.('H3');

Ubase=@(u3)   (s.('U2')+Omega3(u3)*(z-s.('H2'))).*ztop +...
              (s.('U1')+Omega2()*z).*zmid+...
              (s.('U1')+Omega1()*z).*zlow;
       
% plot the velocity profile
ss=get(0,'ScreenSize'); % [left, bottom, width, height]
figure('Position',[ss(3)*0.25,ss(4)*0.25,ss(3)*0.5, ss(4)*0.4])

miny=min(Ubase(s.('U3')));
maxy=max(Ubase(s.('U3')));

plot(Ubase(s.('U3')), z, 'LineWidth',1.5,'Color','red')
hold on
line([miny,maxy],[0,0])
hold on
line([miny,maxy],[s.('H2'),s.('H2')])
hold on
line([miny,maxy],[s.('H2')+s.('H3'),s.('H2')+s.('H3')])
xlim([miny,maxy])
ylim([min(z),max(z)])
hold off
xlabel('velocity')
ylabel('height')
       
% function defining the equation for c       
Sigma1=@(c)     (s.('U0')-c)*(s.('U1')-c)-s.('g12')*s.('H1');
Sigma3=@(u3, c) (s.('U2')-c)*(u3-c)-s.('g32')*s.('H3');
Psi=@(c)        (s.('U2')-c)*(s.('U1')-c);
Fc=@(u3, c)     (s.('H12')*Psi(c)+Sigma1(c))*(s.('H32')*Psi(c)+Sigma3(u3, c))-...
                 s.('H12')*s.('H32')*Psi(c)^2;
cU0=@(u3,c) s.('U0')-c;
cU1=@(u3,c) s.('U1')-c;
cU2=@(u3,c) s.('U2')-c;
cU3=@(u3,c) u3-c;

% plot this function when varying U3
figure('Position',[ss(3)*0.25,ss(4)*0.25,ss(3)*0.5, ss(4)*0.4])
patch([s.('U2'), U3from,  U3from, U3to],...
      [s.('U2'), s.('U2'),U3to,   U3to],lightgray);
hold on
patch([s.('U1'), U3from, U3to,   U3to],...
      [s.('U1'), U3from, U3from, s.('U1'),],lightgray);
hold on

co_fc=fcontour(Fc,[U3from,U3to,U3from,U3to],'LineWidth',1.5,'MeshDensity',MESH_DPS);
set(co_fc,'LevelList',0)
set(co_fc,'LineColor',greybrown)
hold on
co_u3=fcontour(cU3,[U3from,U3to,U3from,U3to],'LineWidth',1.0,'MeshDensity',MESH_DPS);
set(co_u3,'LevelList',0)
set(co_u3,'LineColor',darkgreen)
hold on
co_u1=fcontour(cU1,[U3from,U3to,U3from,U3to],'LineWidth',1.0,'MeshDensity',MESH_DPS);
set(co_u1,'LevelList',0)
set(co_u1,'LineColor','red')
hold on
co_u2=fcontour(cU2,[U3from,U3to,U3from,U3to],'LineWidth',1.0,'MeshDensity',MESH_DPS);
set(co_u2,'LevelList',0)
set(co_u2,'LineColor',orange)
hold on
co_u0=fcontour(cU0,[U3from,U3to,U3from,U3to],'LineWidth',1.0,'MeshDensity',MESH_DPS);
set(co_u0,'LevelList',0)
set(co_u0,'LineColor',blu)
hold on
legend([co_u0,co_u1,co_u2,co_u3],{'$c=U_0$','$c=U_1$','$c=U_2$','$c=U_3$'},'interpreter','latex','Location','NorthEast');
box on

xlabel('$U_3$','interpreter','latex','FontSize',fontsz)
ylabel('$c$','interpreter','latex','FontSize',fontsz)

% edit this part based on where you see a kissing point in the picture!
U3rng_kp=U3rng(-0.7<=U3rng & U3rng<=-0.45);   % area in U3 around where the KP is
ctop=0.2;  % starting value closer to the larger c value
cbot=0.1;  % starting value closer to the smaller c value

c1c2s=zeros(length(U3rng_kp),2);
c1mc2=zeros(length(U3rng_kp),1);
cc=0;   % counter for U3 values for which c1 and c2 to be calculated

for u=U3rng_kp
    cc=cc+1;
    f_tmp=@(c) Fc(u, c);
    c1c2s(cc,1)=fzero(f_tmp,ctop);
    c1c2s(cc,2)=fzero(f_tmp,cbot);
    c1mc2(cc)=abs(c1c2s(cc,1)-c1c2s(cc,2));
end

% find index where is the minimal difference |c1-c2|
idx=find(c1mc2==min(c1mc2));
kp_u3=U3rng_kp(idx);
kp_c1=c1c2s(idx,1);
kp_c2=c1c2s(idx,2);
kp_cbar=(kp_c1+kp_c2)/2;
fprintf(1,'\n*** found a kissing point:\n');
fprintf(1,'c1:%10.4f\nc2:%10.4f\nU3:%10.4f\n',kp_c1,kp_c2,kp_u3);

%plot(kp_u3,kp_c1,'Marker','o','MarkerFaceColor','red','MarkerEdgeColor','red')
%hold on
%plot(kp_u3,kp_c2,'Marker','o','MarkerFaceColor','red','MarkerEdgeColor','red')

text(-0.58,0.22,'$K_3$','Interpreter','latex','FontSize',fontsz)
text(-0.057,-0.43,'$K_4$','Interpreter','latex','FontSize',fontsz)

%%% PART II: EVALUATE INTEGRAL COEFFICIENTS

% modal function and its derivative

b2=1;
b1=@(c)    -b2/(Omega1()*(s.('U0')-c));
m2=@(c)     b2*Sigma1(c)/((s.('U0')-c)*(s.('U1')-c));
m1=@(c)    -b2*s.('H1')/((s.('U0')-c)*(s.('U1')-c))-m2(c)/(Omega2()*(s.('U1')-c));
t2=@(c,u3)  m2(c)*(s.('U2')-c)*(u3-c)/Sigma3(u3, c);
t1=@(c,u3) -t2(c,u3)/(Omega3(u3)*(u3-c));

phi=@(c,u3)...
    (t1(c,u3)+t2(c,u3)./(Omega3(u3)*(s.('U2')-c+Omega3(u3)*(z-s.('H2'))))).*ztop+...
                         (m1(c)+m2(c)./(Omega2()*(s.('U1')-c+Omega2()*z))).*zmid+...
                            (b1(c)+b2./(Omega1()*(s.('U1')-c+Omega1()*z))).*zlow;

phi_z=@(c,u3) -(t2(c,u3)./(s.('U2')-c+Omega3(u3)*(z-s.('H2'))).^2).*ztop-...
                               (m2(c)./(s.('U1')-c+Omega2()*z).^2).*zmid-...
                                  (b2./(s.('U1')-c+Omega1()*z).^2).*zlow;

% check the orthogonality condition; approximate integral using the
% composite trapezoidal rule

kfs=[(1/2) ones(1,nc-2) (1/2)];

% storage for coefficient values

nU3kp=length(U3rng_kp);
npots=1:nU3kp;

storage_c1=zeros(size(npots));
storage_c2=zeros(size(npots));
storage_lambda1=zeros(size(npots));
storage_lambda2=zeros(size(npots));
storage_mu1=zeros(size(npots));
storage_mu2=zeros(size(npots));
storage_sigma=zeros(size(npots));
storage_nu1=zeros(size(npots));
storage_nu2=zeros(size(npots));
storage_gamma1=zeros(size(npots));
storage_gamma2=zeros(size(npots));

for i=1:size(c1c2s,1)
    
    %fprintf(1,'c pair %i: c1=%f c2=%f\n',i,c1c2s(i,1),c1c2s(i,2));
    
    kp_c1=c1c2s(i,1);
    kp_c2=c1c2s(i,2);
    kp_cbar=(kp_c1+kp_c2)/2;
    kp_u3=U3rng_kp(i);
    

    orthc=sum(dz*kfs.*(Ubase(kp_u3)-kp_cbar).*phi_z(kp_c1,kp_u3).*phi_z(kp_c2,kp_u3));
    
    if i==idx
        fprintf(1,'\n*** orthogonality condition:\n%10.7f\n',orthc);
    end
    
    % now calculate coefficients

    l=struct(...
    'alpha1', -2*sum(dz*kfs.*(Ubase(kp_u3)-kp_cbar).*phi_z(kp_c1,kp_u3).^2),...
    'alpha2', -2*sum(dz*kfs.*(Ubase(kp_u3)-kp_cbar).*phi_z(kp_c2,kp_u3).^2),...
    'lambda1',   sum(dz*kfs.*(Ubase(kp_u3)-kp_cbar).^2.*phi(kp_c1,kp_u3).^2),...
    'lambda2',   sum(dz*kfs.*(Ubase(kp_u3)-kp_cbar).^2.*phi(kp_c2,kp_u3).^2),...
    'mu1', (3/2)*sum(dz*kfs.*(Ubase(kp_u3)-kp_cbar).^2.*phi_z(kp_c1,kp_u3).^3),...
    'mu2', (3/2)*sum(dz*kfs.*(Ubase(kp_u3)-kp_cbar).^2.*phi_z(kp_c2,kp_u3).^3),...
    'sigma',     sum(dz*kfs.*(Ubase(kp_u3)-kp_cbar).^2.*phi(kp_c1,kp_u3).*phi(kp_c2,kp_u3)),...
    'nu1', (3/2)*sum(dz*kfs.*(Ubase(kp_u3)-kp_cbar).^2.*phi_z(kp_c1,kp_u3).*phi_z(kp_c2,kp_u3).^2),...
    'nu2', (3/2)*sum(dz*kfs.*(Ubase(kp_u3)-kp_cbar).^2.*phi_z(kp_c1,kp_u3).^2.*phi_z(kp_c2,kp_u3)));
    tmp=(Ubase(kp_u3)-kp_cbar).^2.*phi_z(kp_c1,kp_u3);
    l.('gamma1')=tmp(1);
    tmp=(Ubase(kp_u3)-kp_cbar).^2.*phi_z(kp_c2,kp_u3);
    l.('gamma2')=tmp(1);

    % if alphas are negative, multiply all coefficients by -1
    assert(l.('alpha1')*l.('alpha2')>0, 'alpha1 and alpha2 MUST be of the same sign!')
    if l.('alpha1')<0
        fldn=fieldnames(l);
        for j=1:numel(fldn)
            l.(fldn{j})=-l.(fldn{j});
        end
    end

    % rescale again to make alphas equal to 1
    a=sqrt(l.('alpha1')/l.('alpha2'));
    sc=kp_c1-kp_cbar;
    q=sc*l.('alpha1');

    ll0=struct(...
        'alpha1', 1,...
        'alpha2',  1,...
        'lambda1', l.('lambda1')/q,...
        'lambda2', (a^2/q)*l.('lambda2'),...
        'sigma',   (a/q)*l.('sigma'),...
        'mu1',  l.('mu1')/q,...
        'mu2',  (a^3/q)*l.('mu2'),...
        'nu1', (a^2/q)*l.('nu1'),...
        'nu2', (a/q)*l.('nu2'),...
        'gamma1', l.('gamma1')/q,...
        'gamma2', (a/q)*l.('gamma2'));
    
    w=1/sqrt(abs(ll0.mu1*ll0.mu2));
    
    ll=struct(...
        'alpha1', ll0.alpha1,...
        'alpha2',  ll0.alpha2,...
        'lambda1', ll0.lambda1,...
        'lambda2', ll0.lambda1,...
        'sigma',   ll0.sigma,...
        'mu1',  ll0.('mu1')*w,...
        'mu2',  ll0.('mu2')*w,...
        'nu1', ll0.('nu1')*w,...
        'nu2', ll0.('nu2')*w,...
        'gamma1', ll0.('gamma1')/w,...
        'gamma2', ll0.('gamma2')/w);
    
    storage_c1(i)=kp_c1;
    storage_c2(i)=kp_c2;
    storage_lambda1(i)=ll.('lambda1');
    storage_lambda2(i)=ll.('lambda2');
    storage_mu1(i)=ll.('mu1');
    storage_mu2(i)=ll.('mu2');
    storage_sigma(i)=ll.('sigma');
    storage_nu1(i)=ll.('nu1');
    storage_nu2(i)=ll.('nu2');
    storage_gamma1(i)=ll.('gamma1');
    storage_gamma2(i)=ll.('gamma2');
end

ss=get(0,'ScreenSize'); % [left, bottom, width, height]
figure('Position',[ss(3)*0.25,ss(4)*0.25,ss(3)*0.5, ss(4)*0.5])

subplot(2,3,1)

klr=[73 68 163]./ 255;

plot(npots, storage_c1,'LineStyle','-','Color',klr)
hold on
plot(npots, storage_c2,'LineStyle','--','Color',klr)
xlim([1,nU3kp])
hold on
plot(idx, storage_c1(idx),'Marker','o','MarkerFaceColor',klr,'MarkerEdgeColor',klr) 
hold on
plot(idx, storage_c2(idx),'Marker','o','MarkerFaceColor',klr,'MarkerEdgeColor',klr) 
xlim([1,nU3kp])
set(gca, 'XTickLabel', [])
h=legend('$c_1$','$c_2$','Location','SouthOutside','Orientation','horizontal');
set(h,'Interpreter','latex')
set(h,'FontSize',12)
legend('boxoff')

subplot(2,3,2)

plot(npots, storage_lambda1,'r',...
     npots, storage_lambda2,'--r')
xlim([1,nU3kp])
hold on
plot(idx, storage_lambda1(idx),'Marker','o','MarkerFaceColor','red','MarkerEdgeColor','red') 
hold on
plot(idx, storage_lambda2(idx),'Marker','o','MarkerFaceColor','red','MarkerEdgeColor','red') 
xlim([1,nU3kp])
set(gca, 'XTickLabel', [])
h=legend('$\lambda_1$','$\lambda_2$','Location','SouthOutside','Orientation','horizontal');
set(h,'Interpreter','latex')
legend('boxoff')
set(h,'FontSize',12)
 
subplot(2,3,3)
plot(npots, storage_mu1,'b',...
     npots, storage_mu2,'--b')
xlim([1,nU3kp])
hold on
plot(idx, storage_mu1(idx),'Marker','o','MarkerFaceColor','blue','MarkerEdgeColor','blue') 
hold on
plot(idx, storage_mu2(idx),'Marker','o','MarkerFaceColor','blue','MarkerEdgeColor','blue') 
xlim([1,nU3kp])
set(gca, 'XTickLabel', [])
h=legend('$\mu_1$','$\mu_2$','Location','SouthOutside','Orientation','horizontal');
set(h,'Interpreter','latex')
legend('boxoff')
set(h,'FontSize',12)
 
 
subplot(2,3,4)
plot(npots, storage_sigma,'m')
xlim([1,nU3kp])
hold on
plot(idx, storage_sigma(idx),'Marker','o','MarkerFaceColor','magenta','MarkerEdgeColor','magenta') 
hold on
% plot(idx, storage_sigma(idx),'Marker','o','MarkerFaceColor','magenta','MarkerEdgeColor','magenta') 
xlim([1,nU3kp])
set(gca, 'XTickLabel', [])
h=legend('$\sigma$','Location','SouthOutside','Orientation','horizontal');
set(h,'Interpreter','latex')
legend('boxoff')
set(h,'FontSize',12)

subplot(2,3,5)

klr=[46 174 139]./ 255;

plot(npots, storage_nu1,'LineStyle','-','Color',klr)
hold on
plot(npots, storage_nu2,'LineStyle','--','Color',klr)
xlim([1,nU3kp])
hold on
plot(idx, storage_nu1(idx),'Marker','o','MarkerFaceColor',klr,'MarkerEdgeColor',klr) 
hold on
plot(idx, storage_nu2(idx),'Marker','o','MarkerFaceColor',klr,'MarkerEdgeColor',klr) 
xlim([1,nU3kp])
set(gca, 'XTickLabel', [])
h=legend('$\nu_{1}$','$\nu_{2}$','Location','SouthOutside','Orientation','horizontal');
set(h,'Interpreter','latex')
legend('boxoff')
set(h,'FontSize',12)

subplot(2,3,6)

klr = [201 18 148]./ 255;

plot(npots, storage_gamma1,'LineStyle','-','Color',klr)
hold on
plot(npots, storage_gamma2,'LineStyle','--','Color',klr)
xlim([1,nU3kp])
hold on
plot(idx, storage_gamma1(idx),'Marker','o','MarkerFaceColor',klr,'MarkerEdgeColor',klr) 
hold on
plot(idx, storage_gamma2(idx),'Marker','o','MarkerFaceColor',klr,'MarkerEdgeColor',klr) 
xlim([1,nU3kp])
set(gca, 'XTickLabel', [])
h=legend('$\gamma_1$','$\gamma_2$','Location','SouthOutside','Orientation','horizontal');
set(h,'Interpreter','latex')
set(h,'FontSize',12)
legend('boxoff')