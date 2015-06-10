function S=bedsheet(tfyr,Nx,Nt,earth)
%BEDSHEET  Solves the coupled problem
%     dH/dt = M + Div ( Gam H^{n+2} |grad h|^{n-1} grad h )
% and
%     2 eta |grad| duV/dt + rho g uV + D grad^4 uV = - rhoi g H
% and
%     uE = int int GE(|r-r'|) rhoi H dx dy
% and
%     h = H + uE + uV
% for an isothermal ice sheet on a deforming bed.  Here H(x,y,t) is the ice
% sheet thickness, h(x,y,t) is the ice sheet surface elevation, uV(x,y,t)
% is the vertical displacement of the bed computed from a flat elastic
% plate over a viscous half-space model, and uE(x,y,t) is the vertical
% displacement of the bed computed from a spherical elastic earth model
% whose Green's function is GE.  The equations for uV, uE, and GE are
% described in
%     Bueler, Lingle, and Brown (2006) "Fast computation of a viscoelastic 
%        deformable Earth model for ice sheet simulations", to appear
%        Ann. Glaciol. 46
% The source of this two-part (uV & uE) viscoelastic earth model is
%     Lingle & Clark (1985) "A numerical model of interactions between a
%     marine ice sheet and the solid earth: Application to a West Antarctic
%     ice stream," J. Geophysical Research 90 (C1) 1100--1114.
% The accumulation history is given by the lambda=5 and f=910/3300
% similarity solution in
%     Bueler et al. (2005) "Exact solutions and numerical verification for
%     isothermal ice sheets," J. Glaciol. 51 no. 173, 291--306.
% (In fact the accumulation stops at t0=40034 years.)  This reference also
% describes the numerical method for ice flow model.  
%
% S=BEDSHEET(tfyr,Nx,Nt,'LC')  Uses Lingle & Clark earth deformation
%     model above.  Calls GEFORCONV and VISCDISC for elastic and tweak,
%     respectively.  Returns structure with fields S.h, S.H, S.bed, S.uE, 
%     S.uV, S.xx, S.yy, S.volume, S.Hmax, S.bedmin (all at final time).
%
% S=BEDSHEET(tfyr,Nx,Nt)  Same as above.
%
% Can also compute the simple isostasy model and the ELRA deforming earth 
% model for ice sheet modeling; see, for example,
%     Greve (2001) "Glacial isostasy: Models for the response of the Earth
%     to varying ice loads," in "Continuum Mechanics and Applications in
%     Geophysics and the Environment," Springer, pp. 307--325.
%
% S=BEDSHEET(tfyr,Nx,Nt,'simple')  Simple isostasy with   bed = -f H
%     and   h = (1-f) H   where f=910/3300.  Shows contour map of error in
%     this case.  S has added fields S.err, S.maxerr, S.centererr.
%
% S=BEDSHEET(tfyr,Nx,Nt,'standard')  ELRA model
%
% Example:
%     tic, I = geforconv(32,32,1500,1500); toc  % several minutes
%     save I_1500_32 I
%     S=bedsheet(60000,32,2000,'simple');  % 10 seconds
%     -S.bedmin/S.Hmax  % gives 0.27576
%     SLC=bedsheet(60000,32,2000,'LC');  % about a minute
%     -SLC.bedmin/SLC.Hmax  % gives 0.34248
%
% See also FASTEARTH, GEFORCONV
% ELB 1/15/06; 10/19/06

if nargin<4, earth='LC'; end
if strcmp(earth,'simple'), etype=0;
elseif strcmp(earth,'standard'), etype=1;
elseif strcmp(earth,'LC'), etype=2;
else, error('earthtype not recognized'), end

% physical constants
n=3;            % Glen constant
SperA=31556926; % seconds per year (i.e. 365.2422 days)
A=1e-16/SperA;  % =3.17e-24  1/(Pa^3 s); (EISMINT value) flow law parameter
rhoi=910;       % kg/m^3; density of ice
g=9.81;         % m/s^2; gravity
rhor=3300;      % kg/m^3; density of mantle
Gam=2*(rhoi*g)^n*A/(n+2); % overall constant in deformation discharge q_f
f=rhoi/rhor;    % isostatic parameter
% constants in similarity soln
lam=5;          % for Test C
alf=(2-(n+1)*lam)/(5*n+3);
bet=(1+(2*n+1)*lam)/(5*n+3);
% form of similarity solution at t0
H0=3600;
R0km=750;
R0=R0km*1000;
% t0=40034 years; see equation (9) in Bueler et al. (2005)
t0 = (bet/Gam) * ( (2*n+1)/((1-f)*(n+1)) )^n * (R0^(n+1)/H0^(2*n+1));
tf=tfyr*SperA;
Lkm=1500;
L=Lkm*1e3;     % computational domain for ice flow: (x,y) in [-L,L] x [L,L]

% setup ice flow numerics
dx=2*L/Nx; dy=dx;
% ice flow grid is (Nx+1) x (Nx+1) points
[xx,yy]=meshgrid(linspace(-L,L,Nx+1),linspace(-L,L,Nx+1));
rr=sqrt(xx.^2+yy.^2);
H=zeros(size(rr));  h=H;
uE=H;  uV=uE;  bed=uE+uV;  % initial condition
dt=tf/Nt;
disp(['  [dx = dy = ' num2str(dx/1000) ' (km);  dt = ' num2str(dt/SperA)...
    ' (a)]'])
% used by divflux
Rx=(dt/(dx)^2); Ry=(dt/(dy)^2);
fx=4*dx; fy=4*dy; n2=n+2; nm=(n-1)/2;
in=2:Nx;
ix=in; iy=in; px=ix+1; mx=ix-1; py=iy+1; my=iy-1;

% setup for non-simple earth models
if etype>0
    mindte=200*SperA;      % steps of at least 200 years
    D=5.0e24;              % N m; flexural rigidity of lithosphere plate
    eta=1.0e21;            % Ps s; viscosity of mantle
    Z=2; Lu=Z*L; Nu=Z*Nx;  % computational domain for earth model
    xu=-Lu+dx:dx:Lu;       % domain is periodic: Nu x Nu pts
    [xxu,yyu]=meshgrid(xu,xu);  rru=sqrt(xxu.^2+yyu.^2);
    % Fourier coefficients of powers of Laplacian
    cx=(pi/Lu)*[0:Nu/2 Nu/2-1:-1:1];
    [ccxx,ccyy]=meshgrid(cx,cx);
    cclap=ccxx.^2+ccyy.^2; % lap coeffs
    cchalf=sqrt(cclap);    % lap^{1/2} coeffs
    ccbih=cclap.^2;        % biharmonic coeffs
    sszz=zeros(Nu,Nu);  Hprev=H;
    sh=(Nu/2)-(Nu/(2*Z))+1:(Nu/2)+(Nu/(2*Z))+1;  % for extracting physical
end

if etype==1
    tau=3000*SperA;
    usn=zeros(Nu,Nu);
elseif etype==2
    part1=2*eta*cchalf;
    uun=zeros(Nu,Nu);
    filename=['I_1500_' num2str(Nx) '.mat'];
    doGE=(exist(filename)==2);
    if doGE
        disp(['  [elastic load response matrix FOUND in ' filename ']'])
        S=load(filename);
    else, disp('  [elastic load response matrix NOT found]'), end
end

% time-stepping loop
tic
tprev=0;
for l=0:Nt-1
    % do isothermal ice sheet simulation
    t=dt*l;  tstar=t+0.5*dt;
    Hn=zeros(size(H));
    M=getaccum(xx,yy,tstar);
    Hn(in,in)=H(in,in) - divflux() + dt*M(in,in);
    H=max(0,Hn); % apply boundary condition
    % occasional check for disaster
    if rem(l,30)==0
        if max(abs(H(floor(Nx/2),:))) > H0*2
            error(['instability (blowup) of ice flow at step ' int2str(l)])
        end
    end
    if etype==0 % simple
        h=(1-f)*H;
        bed=-f*H;
    else
        dte=t+dt-tprev;  % time step for earth model
        if dte>mindte  % then take a earth deformation step
            sszz(sh,sh)=-rhoi*g*(H+Hprev)/2;
            if etype==1 % standard
                w=real(ifft2( fft2(sszz)./(rhor*g + D*ccbih) ));
                ra=dte/(2*tau); 
                usn=( (1-ra)*usn + (dte/tau)*w )/(1+ra);
                bed=usn(sh,sh);
            elseif etype==2 % LC
                % for uV:
                part2=(dte/2)*( rhor*g + D*ccbih ); % coeffs in PDE
                right = part1 - part2;  % operator on right of eqn
                left = part1 + part2;   % ... on left of eqn
                % apply Fourier spectral method:
                frhs=right.*fft2(uun) + fft2(dte*sszz);
                uun=real(ifft2( frhs./left ));
                % tweak: see fastearth.m; assume disc radius 1000km
                uun=uun-( sum(uun(1,:))+sum(uun(:,1)) )/(2*Nu);
                H0disc=dx*dx*sum(sum(H))/(pi*1e6^2);  % trapezoid rule
                uun=uun+viscdisc(H0disc,1000,t+dt,Lu/1000);
                uV=uun(sh,sh);
                % for uE: use (convolution) elastic LRM if present
                if doGE,  uE=rhoi*conv2(H,S.I,'same');  end
                bed=uE+uV;
            end

            tprev=t+dt;
            Hprev=H;
        end
        h=H+bed;
    end
end;
time=toc;
disp(['actual comp. time = ' num2str(time) ' secs  (' num2str(time/60) ...
    ' minutes, ' num2str(time/3600) ' hours)'])

% return values; display some
S.h=h; S.H=H; S.bed=bed; S.uE=uE; S.uV=uV; S.xx=xx; S.yy=yy;
S.time=time;
S.Hmax=max(max(H));
disp(['maximum H (at tf)             = ' num2str(S.Hmax) ' (m)']);
S.bedmin=min(min(bed));
disp(['minimum bed (at tf)           = ' num2str(S.bedmin) ' (m)']);
S.ratio=-S.bedmin/S.Hmax;
S.lastRR=Rx*Gam*max(max( Hbr.^n2.*dhsr.^nm ));
disp(['final stability index (at tf) = ' num2str(S.lastRR)]);
S.volume=dx*dy*sum(sum(H));
disp(['final (numerical) ice volume  = ' num2str(S.volume/(1e9*1e6),6) ...
    ' (10^6 km^3)'])

% plot final state
figure(1), clf, subplot(2,1,1)
surf(xx/1000,yy/1000,H);
Hmax=4000;
axis([-Lkm Lkm -Lkm Lkm 0 Hmax]); view(90,0)
xlabel('x (km)'); ylabel('y (km)'); zlabel('H (m)');
title(['Numerical thickness at t=' num2str(tfyr) ' yrs.'])
subplot(2,1,2)
mesh(xx/1000,yy/1000,h);
hold on, mesh(xx/1000,yy/1000,bed); hold off
axis([-Lkm Lkm -Lkm Lkm 1.3*min(min(bed)) 1.2*max(max(h))]); view(90,0)
xlabel('x (km)'); ylabel('y (km)'); zlabel('h (m)');
title(['Numerical surface and bed elevation at t=' num2str(tfyr) ' yrs.'])

% contour of error at final time if 'simple'
if etype==0
    figure(2), clf
    [M,Hexact]=getaccum(xx,yy,tf);  S.err=abs(H-Hexact);
    [Cont,hand] = contour(xx/1000,yy/1000,S.err,[1 5 20 100 500]);
    hand=clabel(Cont,hand);
    axis square, xlabel('x (km)');  ylabel('y (km)');
    title('Absolute thickness error (m).');
    S.maxerr=max(max(S.err));
    disp(['max thickness err (at tf)     = ' num2str(S.maxerr) ' (m)']);
    mid=floor(Nx/2)+1;  S.centererr=S.err(mid,mid);
    disp(['center thickness err (at tf)  = ' num2str(S.centererr) ' (m)']);
end

%%%%% nested subfunctions
    function [M,H]=getaccum(x,y,t)
        % Compute H from Test C, but with simple isostacy added, in
        % Bueler et al (2005) "Exact solutions ... isothermal ice sheets".
        % Stop growth at t0=40033 years.  Multiply by lam t^{-1} to get
        % accumulation; if after t0 accumulation is zero.  After t0
        % solution is essentially Halfar.
        if t<=t0
            ts=t/t0;
            rscl=(ts^(-2))*(sqrt(x.^2+y.^2)/R0);
            H=H0*ts*max(0, 1-rscl.^(4/3) ).^(3/7);
            M=lam*H/t;
        else
            t0post=(t0/2)*(1/18); t=t-t0+t0post; % reset to Halfar w. f
            ts=t/t0post;
            rscl=(ts^(-1/18))*(sqrt(x.^2+y.^2)/R0);
            H=H0*(ts^(-1/9))*max(0, 1-rscl.^(4/3) ).^(3/7);
            M=zeros(size(H));
        end
    end % function getaccum

    function dQ=divflux()
        Hin=H(ix,iy);
        % staggered grid values for H
        Hbr=(H(px,iy) + Hin)/2; Hbl=(Hin + H(mx,iy))/2;
        Hbu=(H(ix,py) + Hin)/2; Hbd=(Hin + H(ix,my))/2;
        % surface gradient, etc., using Mahaffy choice
        hin=h(ix,iy);
        dhsr=((h(px,iy)-hin)/dx).^2 + ((h(px,py)+h(ix,py)-h(px,my)-h(ix,my))/fy).^2;
        dhsl=((hin-h(mx,iy))/dx).^2 + ((h(mx,py)+h(ix,py)-h(mx,my)-h(ix,my))/fy).^2;
        dhsu=((h(px,py)+h(px,iy)-h(mx,py)-h(mx,iy))/fx).^2 + ((h(ix,py)-hin)/dy).^2;
        dhsd=((h(px,iy)+h(px,my)-h(mx,iy)-h(mx,my))/fx).^2 + ((hin-h(ix,my))/dy).^2;
        % form divergence of flux and then mult by dt
        dQx=(Hbr.^n2.*dhsr.^nm).*(h(px,iy)-hin) - (Hbl.^n2.*dhsl.^nm).*(hin-h(mx,iy));
        dQy=(Hbu.^n2.*dhsu.^nm).*(h(ix,py)-hin) - (Hbd.^n2.*dhsd.^nm).*(hin-h(ix,my));
        dQ=-Gam*(Rx*dQx + Ry*dQy);
    end % function divflux

end % function bedsheet
