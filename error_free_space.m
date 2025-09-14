%frequency and wavelenght grid
f0 = 10^10;
GridDiv = 10;
GridXNormal=100;
GridX = 201;
GridY = GridX;
T0=1/f0;

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

EzFreeSpace = zeros(1,GridX-GridXNormal,ceil((12*T0)/dt)+1) ;
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

 
for t = 0:dt:12*T0
 
    Ez(2:end-1,2:end-1)= ((me(2:end-1, 2:end-1) -sigma(2:end-1, 2:end-1).*(dt/2))./(me(2:end-1, 2:end-1) +sigma(2:end-1, 2:end-1).*(dt/2))).*Ez(2:end-1,2:end-1) +...
        (dt./(me(2:end-1 ,2:end-1)+sigma(2:end-1, 2:end-1).*(dt/2))).*((Hy(2:end,2:end-1) - Hy(1:end-1,2:end-1))./dx - (Hx(2:end-1,2:end) - Hx(2:end-1,1:end-1))./dy);
    
   
    if (t>=0) && (t<=10*T0)
    Ez(100,100)=sin(w*t);
    end
    
    EzFreeSpace(1,:,runTimes)=Ez(51,51:151);
  % Update magnetic fields
    Hx(2:end-1,:) = Hx(2:end-1,:) - (dt/(u0*dy)).*(Ez(2:end-1,2:end) - Ez(2:end-1,1:end-1));
    Hy(:,2:end-1) = Hy(:,2:end-1) + (dt/(u0*dx)).*(Ez(2:end,2:end-1) - Ez(1:end-1,2:end-1));

    %some magic to make this work with boundaries


    % Update electric field
    
    
    %temp=Ez(50,50);
    
    

%     Hx
%     Hy
%     Ez


    % Add sinusoidal source to electric field at source location
    

    
    
    
%     Visualize fields at every time step
%     %imagesc(x,y,Ez');
%     surf(x,y,Ez);
%     xlabel('x (m)');
%     ylabel('y (m)');
%     colorbar;
%     title(sprintf('Time: %0.1e s',t));
%     pause(0.0001)
%     runTimes=runTimes+1
end



















