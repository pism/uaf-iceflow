function Hout=verifEIS(type, Tkyr, tfyr, Nx, Mt);
%VERIFEIS  Computes numerical solution of shallow ice equation
%       H_t = M - Div q_f
%    for EISMINT experiments.  Numerical method: type I 
%    explicit finite difference.  Takes advantage of 4-fold symmetry..
%
%verifEIS(type, Tkyr, tf, Nx, Mt);
%   type   = 'fix', 'move' for fixed or moving margin experiments
%   Tkyr   = -1 for no oscillation; 20 or 40 are EISMINT sinusoidal cases
%   tf     = final time (years; run is from t=0 to t=tf)
%   Nx     = number of grid intervals in both x and y
%   Mt     = number of time steps (dt=(tf-t0)/Mt=tadv/Mt)
%
%Notes: 
%   (1) Displays in figure 1.
%   (2) References: (I) Huybrechts, et al. (1996), "The EISMINT benchmarks 
%       testing ice-sheet models," Ann. Glaciol. 23
%       (II) Bueler et al. (2004), "Exact solutions and the verification of 
%       numerical models for isothermal ice sheets", preprint.
%
%Examples:  
%   >> verifEIS('fix',-1,50000,15,5000);    % 8 secs  (EISMINT fixed margin)
%   >> verifEIS('fix',-1,50000,30,20000);   % 63 secs
%   >> verifEIS('fix',-1,50000,60,80000);   % 8.8 mins
%   >> verifEIS('fix',-1,50000,120,320000); % 2.1 hours
%   >> verifEIS('fix',20,200000,15,20000);  % 50 secs  (EISMINT 20ka sinu)
%   >> verifEIS('fix',40,200000,15,20000);  % 47 secs  (EISMINT 40ka sinu)
%   >> verifEIS('move',-1,50000,15,4000);   % 3 secs  (EISMINT moving margin)
%(ELB 4/24/04)

clear H
global H dx dy fx fy n2 nm Rx Ry 
global rho g Gam L n M0 Cs 
global etamonut bsflag r1 r2 theta1 theta2 mustgx mustgy

% physical constants, etc
SperA=31556926; % seconds per year (i.e. 365.2422 days)
A=1e-16/SperA;  %=3.17e-24  1/(Pa^3 s); (EISMINT value) flow law parameter
rho=910; % kg/m^3; density of ice
g=9.81; % m/s^2; gravity
Gam=2*(rho*g)^n*A/(n+2); % overall constant in deformation discharge q_f
tf=tfyr*SperA; TT=Tkyr*1000*SperA;
n=3;
L=750000; 
errcontours=[-500 -200 -100 -70 -50 -30 -20 -10 -5 -1 ...
             1 5 10 20 30 50 70 100 150 500];
Hblowup=10000; % if H reaches this then assume instability

% start numeric comparison
dx=L/Nx; dy=dx; dt=tf/Mt; Rx=(dt/(dx)^2); Ry=(dt/(dy)^2);
[xx,yy]=ndgrid(linspace(0,L,Nx+1),linspace(0,L,Nx+1)); % grid in space
% ndgrid makes coord sys left-handed; better for computation
if false
% Note: Despite p. 4 of EISMINT paper, the initial condition on
%    the sinusoidal runs *is not* the steady result.  Thus there is no need to
%    call recursively to get initial condition.  This is the mechanism.
%if TT>0
   disp('CALLING verifEIS TO GET INITIAL CONDITION (steady)')
   H=verifEIS(type,-1,25000,Nx,ceil(2500*(Nx/15)^2));
   disp('NOW STARTING REQUESTED RUN')
else
   H=zeros(size(xx)); % initial condition is zero
end
t=linspace(0,tf,Mt+1); 
in=2:Nx; fx=4*dx; fy=4*dy; n2=n+2; nm=(n-1)/2;
outice=(xx==L)|(yy==L); % outside of ice if true
disp(['dx   =   dy       = ' num2str(dx/1000) ' km'])
disp(['dt                = ' num2str(dt/SperA) ' years'])

%accumulation setup
if strcmp(type,'fix'),  M0=(0.3/SperA)*ones(size(H)); M=M0;
elseif strcmp(type,'move')
   rr=sqrt(xx.^2+yy.^2);  Rel=450000;
   M=min((0.5/SperA)*ones(size(H)),(10^(-2)/(1000*SperA))*(Rel-rr));
else, error('type must be "fix" or "move"'), end
% figure, mesh(xx,yy,M), axis([0 L 0 L -5/SperA 1/SperA]), return

% time-stepping loop
wbhandle=waitbar(0,'COMPUTING NUMERICAL APPROXIMATION.  Ctrl-C halts.'); tic
for l=1:Mt
   Hn=zeros(size(H));
   % H(1,:) is edge with x=0; H(:,1) is edge with y=0; H(1,1) is corner (x,y)=(0,0)
   
   if TT>0
      t=(l-1)*dt;
      if strcmp(type,'fix')
         M = M0 + (0.2/SperA)*sin(2*pi*t/TT)*ones(size(H));
      else
         Rel=(450+100*sin(2*pi*t/TT))*1000;
         M = min((0.5/SperA)*ones(size(H)),(10^(-2)/(1000*SperA))*(Rel-rr));
      end
   end
   Hn(in,in)=H(in,in) + dt*M(in,in) - divQf(in,in);
   Hn(1,in)=H(1,in) + dt*M(1,in) - divQf(1,in);
   Hn(in,1)=H(in,1) + dt*M(in,1) - divQf(in,1);   
   Hn(1,1)=H(1,1) + dt*M(1,1) - divQf(1,1);

   % apply boundary condition
   if strcmp(type,'fix'), Hn(outice)=0; 
   else, Hn=max(0,Hn); Hn(outice)=0; end
   H=Hn;

   % stability diagnostic (uses only interior points; may miss real max by a bit)
   if l==Mt
      disp(['final stability   = ' num2str(stabindex(in,in)) ' (= max(D)*dt/dx^2)']), end
   % try to check for disaster; waitbar
   if rem(l,30)==0
      if max(abs(H(floor(Nx/2),:))) > Hblowup
         close(wbhandle), error(['instability (blowup) at step ' int2str(l)]), end
      waitbar(l/(Mt+1)), end
   % if lots of steps, estimate compute time
   if rem(l,1500)==0
      remsecs=ceil( (Mt+1-l) * (toc/l) ); close(wbhandle)
      wbhandle=waitbar(0, ['COMPUTING ... ESTIMATED WAIT TIME ' ...
            int2str(remsecs) ' secs']);
      waitbar(l/(Mt+1),wbhandle), end
end;
disp(['actual comp. time = ' num2str(toc) ' secs'])
close(wbhandle)

% report: plot numerical final state; give values
figure(1), clf
surf(xx/1000,yy/1000,H);
if strcmp(type,'fix'), axis([0 L/1000 0 L/1000 0 3800]); 
else axis([0 L/1000 0 L/1000 0 3200]); end
view(90,0), xlabel('x in km'); ylabel('y in km'); zlabel('h in m');
title(['Numerical final state at t = ' num2str(tf/SperA) ' yrs.  Rotatable 3D fig.'])
disp(['central height    = ' num2str(H(1,1))]);
mm=floor(400000/dx+10*eps);  midflux=((Qfx(mm,1)+Qfx(mm+1,1))/2);
disp(['midpoint flux     = ' num2str(midflux*SperA/100) '    10^2 m^2/a']);
if strcmp(type,'move')
   pos=max(max(xx(H>0)))+dx;
   disp(['pos of margin     = ' num2str(pos/1000) ' km'])
end
Hout=H;

%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%
function dQ=divQf(ix,iy)
% (numerical) Divergence of deformation (flow) flux
global H dx dy fx fy n2 nm Rx Ry Gam
if length(ix)==1, px=2; mx=2; else, px=ix+1; mx=ix-1; end
if length(iy)==1, py=2; my=2; else, py=iy+1; my=iy-1; end
Hin=H(ix,iy);
Hbr=(H(px,iy) + Hin)/2; Hbl=(Hin + H(mx,iy))/2;
Hbu=(H(ix,py) + Hin)/2; Hbd=(Hin + H(ix,my))/2;
dHsr=((H(px,iy)-Hin)/dx).^2 + ((H(px,py)+H(ix,py)-H(px,my)-H(ix,my))/fy).^2;
dHsl=((Hin-H(mx,iy))/dx).^2 + ((H(mx,py)+H(ix,py)-H(mx,my)-H(ix,my))/fy).^2;
dHsu=((H(px,py)+H(px,iy)-H(mx,py)-H(mx,iy))/fx).^2 + ((H(ix,py)-Hin)/dy).^2;
dHsd=((H(px,iy)+H(px,my)-H(mx,iy)-H(mx,my))/fx).^2 + ((Hin-H(ix,my))/dy).^2;
dQx=(Hbr.^n2.*dHsr.^nm).*(H(px,iy)-Hin) - (Hbl.^n2.*dHsl.^nm).*(Hin-H(mx,iy));
dQy=(Hbu.^n2.*dHsu.^nm).*(H(ix,py)-Hin) - (Hbd.^n2.*dHsd.^nm).*(Hin-H(ix,my));
dQ=-Gam*(Rx*dQx + Ry*dQy);

function Q=Qfx(ix,iy)
% (numerical) deformation flux in x dir; only used for diagnostic at end
global H dx dy fx fy n2 nm Rx Ry Gam
if ix==1, px=2; mx=2; else, px=ix+1; mx=ix-1; end
if iy==1, py=2; my=2; else, py=iy+1; my=iy-1; end
Hin=H(ix,iy);
Hbr=(H(px,iy) + Hin)/2;
dHsr=((H(px,iy)-Hin)/dx).^2 + ((H(px,py)+H(ix,py)-H(px,my)-H(ix,my))/fy).^2;
Q=-Gam*(Hbr.^n2.*dHsr.^nm).*(H(px,iy)-Hin)/dx;

function RR=stabindex(ix,iy)
% compute certain diffusivities for max estimate
global H dx dy fx fy n2 nm Rx Ry Gam rho g bsflag mustgx
px=ix+1; py=iy+1; my=iy-1;
Hin=H(ix,iy); Hbr=(H(px,iy) + Hin)/2;
dHsr=((H(px,iy)-Hin)/dx).^2 + ((H(px,py)+H(ix,py)-H(px,my)-H(ix,my))/fy).^2;
Df=Gam*max(max( Hbr.^n2.*dHsr.^nm )); RR=Rx*Df;
if bsflag
   Db=rho*g*max(max( mustgx(ix,iy).*Hbr.^2 )); RR=RR+Rx*Db;
end

