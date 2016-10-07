function[C,c,n,h,x,y]=IMEX_Central_Diff()
% Isotrop ohne Reakion (IMEX)
%function[timee]=IMEX_Central_Diff()
% Number of steps and cycles
Steps=50;
r=0;
Zyklen=200;
% Time step
tau=0.02;
% Length
L=0.002;


% Diffusion coefficient
diffusion_0=1*10^-10;

% F�r Geschwindigkeit
%Viskosit�t
nu=10;
% Permeabilit�t

k_0=7.5*10^-16;
% Fl�che in Meter
area=10^-12;
% Konstante M
konstant2=4.3;
% Void ratio = Porosit�t
phi_0=0.8;
%Dehnung �ber Zeit tr(e)
e=linspace(0,0.1,Steps);
% Druck gradient
delta_p=linspace(0,2*10^6,Steps);
delta_p_Abstieg=linspace(2*10^6,0,Steps);
% Anstieg
e_Anstieg=linspace(0.1,0,Steps);

%a=10^-3;
%Setting up the grid

h=L/75;
x = 0:h:L;
y = 0:h:L;
n=length(x);
ndgrid(n,n); 

%concentration matrix
c=zeros(n^2,1);
c_0=1;

Ident=speye(n^2,n^2);
 Ident_n=speye(n,n);
      %% Building matrix for Laplace Operator
    T= (1/h^2)*(gallery('tridiag',n,1,-2,1));   
    A=kron(T,Ident_n)+kron(Ident_n,T);
   
   %% building gradient matrix (backwards scheme) upwind
      A1 = (1/(2*h))*(gallery('tridiag',n^2,-1,0,1)); 
%% Belastung
tic
while(r<Zyklen)  
for i=1:Steps
   
    % porosity at step
    phi=1-(1-phi_0)*(1+e(i));
    % velocity darcy law
    k=k_0*exp(konstant2*(phi-phi_0)/(1-phi));    
    k=k/nu;
    Q=k*area*delta_p(i)/L;
    q=Q/(phi*area);
   
    
    
    % diffusion
    diffusion=diffusion_0*exp(10*(phi-phi_0)/(phi*phi_0));
    
   
      

    

    % M= (I - tau*q*D^-*0.5)
    M=  Ident + tau*q/2*A1  ;
    
    % N=  (I +tau(-a + D^+D^-) + tau/2 * (qD^-)
    N = Ident + tau*diffusion*A  - (tau/2)*A1*q;

  %% Boundary conditions
    rhs=N*c;
    [c]=Boundary(M,rhs,n,c_0,h,q,Q);

    
end
%Static
for i=1:Steps

    

  
  

 % M= (I - tau*q*D^-*0.5)
    M=  Ident + tau*q/2*A1  ;
    
    % N=  (I +tau(-a + D^+D^-) + tau/2 * (qD^-)
    
    N = Ident + tau*diffusion*A  - (tau/2)*A1*q;
  %% Boundary conditions
    rhs=N*c;
    [c]=Boundary(M,rhs,n,c_0,h,q,Q);
    
end
%% Entlastung
for i=1:Steps
    % porosity at step
    phi=1-((1-phi_0)*(1+e_Anstieg(i)));
    % velocity darcy law
     k=k_0*exp(konstant2* (phi-phi_0)/(1-phi));    
     k= k/nu;
     Q=k*area*delta_p_Abstieg(i)/L;
     q=Q/(phi*area);
 
    % diffusion
    diffusion=diffusion_0*exp(10*(phi-phi_0)/(phi*phi_0));
    
   

 % M= (I + tau*q*D^-*0.5)
    M=  Ident + tau*q/2*A1  ;
    
    % N=  (I -tau(-a + D^+D^-) - tau/2 * (qD^-)
    N = Ident + tau*diffusion*A  - (tau/2)*A1*q;
    
  %% Boundary conditions
    rhs=N*c;
    [c]=Boundary(M,rhs,n,c_0,h,q,Q);
 
end
   
    r=r+1;
    r
end
 toc  
   C=reshape(c,[n,n]);
   C=C';
   x=x*1000;
   y=fliplr(y);
   y=y*1000;
   [X,Y] = meshgrid(x,y); 
   pcolor (X,Y,C);
   shading interp;
   colormap(jet(4096));
   title ("Konvektion Diffusion im isotropen Medium","fontsize", 14,"fontweight","bold");
   xlabel("x [mm]","fontsize", 13,"fontweight","bold");
   colorbar();
   ylabel("y [mm]","fontsize", 13,"fontweight","bold");
   set(gca,'ytick',[0 1 2]);
   cbh = findobj( gcf(), 'tag', 'colorbar')
     set( cbh,                %now change some colorbar properties
    'linewidth', 1.2,
    'tickdir', 'out',
    'ticklength',[0.005,0.005],
    'title','Konzentration [-]')
   
   print -depsc IMEX_Conv_Diff_Iso.eps
end



