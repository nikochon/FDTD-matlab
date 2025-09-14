%frequency and wavelenght grid
f0 = 10^10;
T0=1/f0;
GridDiv = 10;
GridX = 101;
GridY = GridX;

%free space constants
e0 = 8.54*10^-12;
u0 = 4*pi*10^-7;
c0 = 3*10^8;
L0 = c0/f0;

Px=zeros(1,GridX);
Py=zeros(1,GridY);


for i= 1:GridX
    Px(1,i)=i;
    
end

for i= 1:GridY
    Py(1,i)=i;
    
end

x=Px.*(L0/GridDiv);
y=Py.*(L0/GridDiv);

%increments

dx=L0/GridDiv;
dy=L0/GridDiv;
dt=5*10^-12;       %anything lower than 7.071*10^-12 to be stable
tmax=dt;        %simulation time


%field initialiatons
Ez = zeros(length(x),length(y)); 
Hx = zeros(length(x),length(y)-1); 
Hy = zeros(length(x)-1,length(y));

mu = zeros(length(x),length(y));
me = zeros(length(x),length(y));
sigma = zeros(length(x),length(y));


%source

fs = 10^10;
w = 2*pi*fs;
temp=0;

runTimes=1; %debug

%cylider
radius=L0;
c_x=50*dx+3*L0;
c_y=50*dy;
c_er=3.2;
c_mu=1;
c_sigma=1.2; %see later


for i=1:GridX
    for j=1:GridY
        if false
            mu(i,j)=c_mu*u0;
            me(i,j)=c_er*e0;
            sigma(i,j)=c_sigma;
        else
            mu(i,j)=u0;
            me(i,j)=e0;
            sigma(i,j)=0;
        end
    end
end

%pml initialiata =====================================
Npml=30;
O=3;
R=10^-6;


sigmaE=zeros(Npml,1);
sigmaHx=zeros(Npml,1);
sigmaHy=zeros(Npml,1);

Hx_pml=zeros(Npml,GridY-1);
Hy_pml=zeros(Npml,GridY);
Ezx_pml=zeros(Npml,GridY);
Ezy_pml=zeros(Npml,GridY);

se= -e0*c0*log(R)/(2^(O+2)*dx*Npml^(O+1));

for i=1:Npml
sigmaE(i)=se*((2*i+1)^(O+1)-(2*i-1)^(O+1)) ;
end

sh=u0/e0*se;

for i=1:Npml
sigmaHx(i)=sh*((2*i+1)^(O+1)-(2*i-1)^(O+1));
sigmaHy(i)=sh*((2*(i+0.5)+1)^(O+1)-(2*(i+0.5)-1)^(O+1));
end

sigmaE=fliplr(sigmaE); 
sigmaHx=fliplr(sigmaHx);
sigmaHy=fliplr(sigmaHy);



Ca_pml=exp(1).^(-sigmaE.*dt/e0);
Cb_pml=(1-Ca_pml)./(sigmaE.*dx);
Dax_pml= exp(1).^( -sigmaHx.*dt./u0);
Day_pml=exp(1).^( -sigmaHy.*dt./u0);
Dbx_pml=(1-Dax_pml)./(sigmaHx.*dx); 
Dby_pml=(1-Day_pml)./(sigmaHy.*dx);

MCa_pml=zeros(length(Ca_pml),GridY);
MCb_pml=zeros(length(Cb_pml),GridY);
MDax_pml=zeros(length(Dax_pml),GridY);
MDbx_pml=zeros(length(Dbx_pml),GridY);
MDay_pml=zeros(length(Day_pml),GridY);
MDby_pml=zeros(length(Dby_pml),GridY);

for j= 1 : GridY

    MCa_pml(:,j)=Ca_pml;
    MCb_pml(:,j)=Cb_pml;
    MDax_pml(:,j)=Dax_pml;
    MDbx_pml(:,j)=Dbx_pml;
    MDay_pml(:,j)=Day_pml;
    MDby_pml(:,j)=Dby_pml;

end


%outer PML Γ=1
Hx_pml(1,:)=0;
Hy_pml(1:end,1)=0;
Hy_pml(1:end,end)=0;

Ezx_pml(1,1:end)=0;
Ezx_pml(1:end,1)=0;
%Ezx_pml(1:end,end)=0;

Ezy_pml(1,1:end)=0;
Ezy_pml(1:end,1)=0;
%Ezy_pml(1:end,end)=0;



EzPML=zeros(1,GridY,ceil((12*T0)/dt)+1);
 
for t = 0:dt:12*T0
 
    
   

    

   
  

 


    % Update electric
    % field========================================================
    
    Ez(2:end-1,2:end-1)= ((me(2:end-1, 2:end-1) -sigma(2:end-1, 2:end-1).*(dt/2))./(me(2:end-1, 2:end-1) +sigma(2:end-1, 2:end-1).*(dt/2))).*Ez(2:end-1,2:end-1) +...
        (dt./(me(2:end-1 ,2:end-1)+sigma(2:end-1, 2:end-1).*(dt/2))).*((Hy(2:end,2:end-1) - Hy(1:end-1,2:end-1))./dx - (Hx(2:end-1,2:end) - Hx(2:end-1,1:end-1))./dy);
    
    Ez(1,2:end-1)=((me(1, 2:end-1) -sigma(1,2:end-1).*(dt/2))./(me(1,2:end-1) +sigma(1, 2:end-1).*(dt/2))).*Ez(1,2:end-1)...
        +(dt./(me(1 ,2:end-1)+sigma(1,2:end-1).*(dt/2))).*((Hy(1,2:end-1)-Hy_pml(end,2:end-1))./dx  - (Hx(1,2:end) - Hx(1,1:end-1))./dy); 



   
    %PML FIELDS========================================================= EZ
   
    
    Ezx_pml(2:end , 2:end-1) = MCa_pml(2:end , 2:end-1).* Ezx_pml(2:end , 2:end-1)...
        +MCb_pml(2:end,2:end-1).*(Hy_pml(2:end,2:end-1)-Hy_pml(1:end-1,2:end-1));

    Ezy_pml(2:end , 2:end-1) = MCa_pml(2:end , 2:end-1).* Ezy_pml(2:end , 2:end-1)...
        +MCb_pml(2:end,2:end-1).*(Hx_pml(2:end,1:end-1)-Hx_pml(2:end,2:end));

    % Add sinusoidal source to electric field at source location
    if (t>=0) && (t<=10*T0)
    Ez(50,50)=sin(w*t);
    end
    EzPML(1,:,runTimes)=Ez(1,:);


    % Update magnetic fields
    Hx(2:end-1,:) = Hx(2:end-1,:) - (dt/(u0*dy)).*(Ez(2:end-1,2:end) - Ez(2:end-1,1:end-1));


    %PML=================================== HX
    Hx_pml(2:end,:) = MDax_pml(2:end,1:end-1).*Hx_pml(2:end,:) ...
        +MDbx_pml(2:end,1:end-1).*( (Ezx_pml(2:end,1:end-1)+Ezy_pml(2:end,1:end-1)) ...
        - (Ezx_pml(2:end,2:end)+Ezy_pml(2:end,2:end) ) );


 %%%%%%%%%%%%%%%%%%%%%%%%%
    Hy(:,2:end-1) = Hy(:,2:end-1) + (dt/(u0*dx)).*(Ez(2:end,2:end-1) - Ez(1:end-1,2:end-1));
    
   
    

        Hy_pml(1:end-1,2:end-1) = MDay_pml(1:end-1,2:end-1).*Hy_pml(1:end-1,2:end-1) ...
        +MDby_pml(1:end-1,2:end-1).*( (Ezx_pml(2:end,2:end-1)+Ezy_pml(2:end,2:end-1)) ...
        - (Ezx_pml(1:end-1,2:end-1)+Ezy_pml(1:end-1,2:end-1) ) );

        Hy_pml(end,:)= MDay_pml(end,:).*Hy_pml(end,:)...
            +MDby_pml(end,:).*(Ez(1,:)-  (Ezx_pml(end,:)+Ezy_pml(end,:) ));

        

    
%     figure(1)
%     surf(x,y,Ez);
%     title('Ez');
%     zlim([-1,1]);
%     xlabel('x/λ');
%     ylabel('y/λ');
%     zlabel('z');

    % hold off;
    drawnow();


    
    

    
    
    
    % Visualize fields at every time step
%     figure(2)
%     imagesc(x,y,Ez');
%     %imagesc(Ezy_pml)
%     %surf(x,y,Ez);
%     xlabel('x (m)');
%     ylabel('y (m)');
%     colorbar;
%     title(sprintf('Time: %0.1e s',t));
%     pause(0.0001)
    runTimes=runTimes+1
end



















