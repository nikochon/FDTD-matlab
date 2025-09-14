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
source_location = [50, 50];
fs = 10^10;
w = 2*pi*fs;
temp=0;

runTimes=0; %debug

%cylider
radius=L0;
c_x=50*dx+3*L0;
c_y=50*dy;
c_er=3.2;
c_mu=1;
c_sigma=1.2; %see later


for i=1:GridX
    for j=1:GridY
        if sqrt((i*dx-c_x)^2 + (j*dy-c_y)^2) <= radius
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

 
for t = 0:dt:10*T0
 
    
   
    % Update electric field
    Ez(2:end-1,2:end-1)= Ez(2:end-1,2:end-1) + (dt/e0).*((Hy(2:end,2:end-1) - Hy(1:end-1,2:end-1))./dx - (Hx(2:end-1,2:end) - Hx(2:end-1,1:end-1))./dy);
    
    
    if (t>=0) && (t<=10*T0)
    Ez(source_location(1),source_location(2))=sin(w*t);
    end
    

   
  % Update magnetic fields
    Hx(2:end-1,:) = Hx(2:end-1,:) - (dt/(u0*dy)).*(Ez(2:end-1,2:end) - Ez(2:end-1,1:end-1));
    Hy(:,2:end-1) = Hy(:,2:end-1) + (dt/(u0*dx)).*(Ez(2:end,2:end-1) - Ez(1:end-1,2:end-1));

    


%     Hx
%     Hy
%     Ez


    
    

    
    
    
    % Visualize fields at every time step
    surf(x,y,Ez'+1)
    hold on
    imagesc(x,y,Ez'+1)
    hold off

    xlabel('x (m)');
    ylabel('y (m)');
    
    title(sprintf('Time: %0.1e s',t));
    pause(0.0001)
    runTimes=runTimes+1
end

    surf(x,y,Ez'+1)
    hold on
    imagesc(x,y,Ez'+1)
    colorbar;



















