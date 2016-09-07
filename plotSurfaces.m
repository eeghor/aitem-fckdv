%%%%
%%%% PLOT SURFACES for CHAPTER 5 of THESIS // 07 Sep, 2016
%%%%

function plotSurfaces(pa,boa,gr,cs,csp,ua,PNT,MIN_z,VIEW_ANGLE)

% all sufraces will be with momentum and parameter b as x and y and
% something else as z

[P,B]=meshgrid(pa,boa);

% how many vlaues of p and b we have

NVAL_p=length(pa);
NVAL_b=length(boa);

% storages

CC=zeros(NVAL_b, NVAL_p);
DC=zeros(NVAL_b, NVAL_p);
GR=zeros(NVAL_b, NVAL_p);
PA=zeros(NVAL_b, NVAL_p);
PB=zeros(NVAL_b, NVAL_p);
Cm1=zeros(NVAL_b, NVAL_p);

cn=0; % counter for points corresponding to c>-1

for p=1:NVAL_p
    for b=1:NVAL_b
        
        if cs(p,b)<=-1   % we are only interested in where c<-1!
            GR(b,p)=real(gr(p,b));
            CC(b,p)=cs(p,b);
            DC(b,p)=csp(p,b);
            Cm1(b,p)=-2; % always -2 for the right c values
            As=squeeze(ua(p,b,1,:));
            Bs=squeeze(ua(p,b,2,:));
            PA(b,p)=sum(As.^2)/p;
            PB(b,p)=sum(Bs.^2)/p;
        else
            GR(b,p)=nan;
            CC(b,p)=nan; 
            DC(b,p)=nan;
            Cm1(b,p)=0; % always 0 for the wrong c values
            cn=cn+1;
            PA(b,p)=nan;
            PB(b,p)=nan;
        end
        
    end
end

% set all values corresponding to momentum equal to 1 to NaNs;
% this is to make sure that lines corresp. to momentum 1 do NOT
% appear in the plots

CC(:,NVAL_p)=nan;
Cm1(:,NVAL_p)=nan;
GR(:,NVAL_p)=nan;
DC(:,NVAL_p)=nan;
PA(:,NVAL_p)=nan;
PB(:,NVAL_p)=nan;

ss=get(0,'ScreenSize'); % [left, bottom, width, height]

figure('Position',[ss(3)*0.2,ss(4)*0.25,ss(3)*0.6, ss(4)*0.4])

subplot(1,2,1)

    % plot the surface b vs momentum vs phase speed c
    
    h=surfc(P,B,CC);
    zlim([MIN_z,-1]);
    xlabel('momentum')
    ylabel('b')
    zlabel('c')
    colormap gray
    brighten(0.6)
    set(h(1),'FaceColor',[0.8 0.8 0.8]) % surface color
    set(h(2),'Visible','off')  % hide contour
    view(VIEW_ANGLE)
    hold on
    
    % add a surface with b vs p vs dc/dp
    
    r=surfc(P,B,DC);
    zlim([MIN_z,-1]);
    set(r(1),'Visible','off')  % hide surface
    set(r(2),'LineWidth',1.5) 
    set(r(2),'LineColor', [22 154 90]./255)  % green; was red/pink [230 11 46]./255
    set(r(2), 'LineStyle','-')
    set(r(2),'LevelListMode','manual')
    set(r(2),'LevelList',0)
    
    % add a surface where c > -1 (if there are such values of c)

    if cn>0
        
        hold on
        j=surfc(P,B,Cm1);  % add a contour where Cm1=0 (corresp. to where c>-1)
        set(j(1),'Visible','off')  % hide surface
        set(j(2),'LineWidth',1.5)
        set(j(2),'LineColor', [108 122 167] ./ 255)
        set(j(2),'LevelListMode','manual')
        set(j(2),'LevelList',0)
        hold on
     
        fprintf(1,'# points corresponding to c > -1: %d\n',cn);
    else
        fprintf(1,'all values of c are <-1\n');
    end
    
    % add a point for which there will an unsteady simulation example
    
    zLimits=get(gca,'ZLim');  % need to put the point on the bottom
    plot3([pa(PNT(1)),pa(PNT(1))],...
          [boa(PNT(2)),boa(PNT(2))],...
          [zLimits(1),cs(PNT(1),PNT(2))],'--ob','MarkerFaceColor','b')

subplot(1,2,2)
    
    % these contour plots considered 3D!
    % first, the dc/dp=0 contour
    contour(DC,[0,0],'LineWidth',1.5,'LineColor',[22 154 90]./255);
    hold on
    
    % now add the growth rate contour
    contour(GR,[1e-9,1e-9],'LineWidth',1.5,'LineColor',[255 3 100]./255);
    
    % and if not all c are <-1, add the contour where c >-1
    if cn>0 
        hold on 
        contour(Cm1,[0,0],'LineWidth',1.5,'LineColor',[108 122 167] ./ 255);
    end
    
    set(gca,'XTick',10:20:NVAL_p)
    set(gca,'XTickLabel',pa(get(gca,'XTick')))
    set(gca,'YTick',1:10:NVAL_b)
    set(gca,'YTickLabel',boa(get(gca,'YTick')))
    xlabel('momentum')
    ylabel('b')
    fn={'dc/dp=0', 'zero growth rate','c>-1'};
    legend(fn,'Location','North');
    grid on
    
    % add the point for which there will be an unsteady simulation example
    
    hold on
    plot3(PNT(1),PNT(2),0,'ob','MarkerFaceColor','b')

% this figure will show the momentum ratios
figure('Position',[ss(3)*0.2,ss(4)*0.25,ss(3)*0.6, ss(4)*0.4])

subplot(1,2,1)
    surf(P,B,2*PA./(PA+PB)-1);
    xlabel('momentum')
    ylabel('b')
    zlabel('(PA-PB)/P')
    colormap gray
    brighten(0.8)
    view(VIEW_ANGLE)

subplot(1,2,2)
    surf(P,B,PB./(PA+PB));
    xlabel('momentum')
    ylabel('b')
    zlabel('PB/P')
    colormap gray
    brighten(0.8)
    view(VIEW_ANGLE)

end

