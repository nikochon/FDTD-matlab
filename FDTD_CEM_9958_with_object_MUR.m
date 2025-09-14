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
dt=5*10^-12;  %anything lower than 7.071*10^-12 to be stable
tTimes=0:dt:12*T0;
      


%field initialiatons
Ez = zeros(length(x),length(y)); 
Hx = zeros(length(x),length(y)-1); 
Hy = zeros(length(x)-1,length(y));
EzMidMUR=(length(tTimes));
EzLeftMUR=(length(tTimes));

mu = zeros(length(x),length(y));
me = zeros(length(x),length(y));
sigma = zeros(length(x),length(y));

    %mur conditions
    Ez1_t1=Ez(1,1:end);
    Ez2_t1=Ez(2,1:end);
    





%source============================
source_location = [50, 50];
fs = 10^10;
w = 2*pi*fs;


runTimes=0; 

%cylider==============================
radius=L0;
c_x=50*dx+3*L0;
c_y=50*dy;
c_er=3.2;
c_mu=1;
c_sigma=1.2; 


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

 

for t = tTimes
 
   
  Ez1_t2=Ez1_t1;  %n-1
  Ez1_t1 = Ez(1,1:end); %n
  Ez2_t2=Ez2_t1; %n-1
  Ez2_t1 = Ez(2,1:end);  %n

   
 % Update electric field=============================================
    Ez(2:end-1,2:end-1)= ((me(2:end-1, 2:end-1) -sigma(2:end-1, 2:end-1).*(dt/2))./(me(2:end-1, 2:end-1) +sigma(2:end-1, 2:end-1).*(dt/2))).*Ez(2:end-1,2:end-1) +...
        (dt./(me(2:end-1 ,2:end-1)+sigma(2:end-1, 2:end-1).*(dt/2))).*((Hy(2:end,2:end-1) - Hy(1:end-1,2:end-1))./dx - (Hx(2:end-1,2:end) - Hx(2:end-1,1:end-1))./dy);
    
    Ez(1,1)=Ez2_t1(1,1)-((dx- c0*dt)/(dx+c0*dt))*(Ez(2,1)-Ez1_t1(1,1));
    Ez(1,end)=Ez2_t1(1,end)-((dx- c0*dt)/(dx+c0*dt))*(Ez(2,end)-Ez1_t1(1,end));



    Ez(1,2:end-1)=-Ez2_t2(1,2:end-1) - ((dx- c0*dt)/(dx+c0*dt))*(Ez(2,2:end-1)+Ez1_t2(1,2:end-1)) ...
       +((2*dx)/(dx+c0*dt))*(Ez1_t1(1,2:end-1)+Ez2_t1(1,2:end-1)) ...
       +((c0^2 * dt^2 * dx)/(2*dx^2*(dx+c0*dt)))...
       *(Ez1_t1(1,3:end)-2*Ez1_t1(1,2:end-1)+Ez1_t1(1,1:end-2)...
       +Ez2_t1(1,3:end)-2*Ez2_t1(1,2:end-1)+Ez2_t1(1,1:end-2));
        
    if (t>=0) && (t<=10*T0)
    Ez(source_location(1),source_location(2))=sin(w*t);
    end
    
  EzMidMUR(runTimes+1)=Ez(50,10);
  EzLeftMUR(runTimes+1)=Ez(10,10);

  % Update magnetic fields==============================================
    Hx(2:end-1,:) = Hx(2:end-1,:) - (dt/(u0*dy)).*(Ez(2:end-1,2:end) - Ez(2:end-1,1:end-1));
    Hy(:,2:end-1) = Hy(:,2:end-1) + (dt/(u0*dx)).*(Ez(2:end,2:end-1) - Ez(1:end-1,2:end-1));

    


   

    
    


%     Hx
%     Hy
%     Ez


    
%     figure(1)
%     surf(x,y,Ez);
%     title('Ez');
%     zlim([-1,1]);
%     xlabel('x/λ');
%     ylabel('y/λ');
%     zlabel('z');
    
    
    
    % Visualize fields at every time step
    surf(x,y,Ez'+1)
    hold on
    imagesc(x,y,Ez'+1) 
    hold off
    
    xlabel('x (m)');
    ylabel('y (m)');
    
    title(sprintf('Time: %0.1e s , Ez',t));

    pause(0.0001)
    runTimes=runTimes+1
end


    surf(x,y,Ez'+1)
    hold on
    imagesc(x,y,Ez'+1)
    colorbar;
















