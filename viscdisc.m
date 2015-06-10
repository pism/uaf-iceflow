function u=viscdisc(H0,R0km,tyrs,rkm)
% VISCDISC  Solution of viscous half-space overlain by rigid lithosphere in
% case of disc load at origin.  Computed by Hankel transform means for 
% comparison in
%     Bueler, Lingle, and Brown (2006) "Fast computation of a viscoelastic 
%        deformable Earth model for ice sheet simulations", to appear
%        Ann. Glaciol. 46
% Also called by FASTEARTH to correct boundary condition.
%
% U=VISCDISC(H0,R0KM,TYRS,RKM)  computes the displacement U at
%    time TYRS (years) at radii RKM (km) caused by a disc of radius
%    R0KM (km) and thickness H0 (m).  RKM may be a vector; U will be a 
%    vector of the same length.
%
% U=VISCDISC(H0,R0KM,'INF',RKM)  computes the equilibrium 
%    displacement.
%
% Example; continues the above.  Plots the Green's function as in figure
% 5(a) in Lingle & Clark (1985), but at fewer times.  Takes a minute or so.
%
%    r=10.^(-2:.2:4);  % 31 r values
%    tt=[20 200 500 2000 5000 20000 50000 100000];
%    R0=0.001;  H0=1/(pi*910);  % small disc with mass 1 kg
%    for j=1:8,  uu(:,j)=viscdisc(H0,0.001,tt(j),r);  end
%    figure,  for j=1:8, semilogx(r,uu(:,j)),  hold on,  end
%    hold off,  xlabel('r (km)'),  ylabel('displacement (m)')
%    axis([r(1) r(end) -3.5e-15 0.5e-15]), grid on
%
% See also FASTEARTH.
% ELB 12/7/05; 10/19/06

g=9.81; % m/s^2
rhoi=910; % kg/m^3
rhor=3300;  % kg/m^3
rg=rhor*g; % combined constant
D=5.0e24; % N m; flexural rigidity of lithosphere
eta=1.0e21;  % Pa s; viscosity of mantle
spera=3.1556926e7;  % seconds per year

R0=R0km*1000;
TOL=eps;
pts=[10.^(-3:-0.05:-10) 1.0e-14];

if (isstr(tyrs)&&strcmp(tyrs,'inf'))
    for k=1:length(rkm)
        rk=rkm(k)*1000;
        result=quadl(@equilgrand,pts(1),100.0*pts(1),TOL); %kap->infty tail
        for j=1:length(pts)-1
            result=result+quadl(@equilgrand,pts(j+1),pts(j),TOL);
        end
        u(k)=-rhoi*g*H0*R0*result;
    end
else
    t=tyrs*spera;
    for k=1:length(rkm)
        rk=rkm(k)*1000;
        result=quadl(@integrand,pts(1),100.0*pts(1),TOL); % kap->infty tail
        for j=1:length(pts)-1
            result=result+quadl(@integrand,pts(j+1),pts(j),TOL);
        end
        u(k)=rhoi*g*H0*R0*result;
    end
end

    function y=integrand(kap)
    % integrand of inverse Hankel transform
    beta=rg + D*kap.^4;
    expdiff=exp(-beta*t./(2*eta*kap))-ones(size(kap));
    y=expdiff.*besselj(1.0,kap*R0).*besselj(0.0,kap*rk)./beta;
    end

    function y=equilgrand(kap)
    % integrand of inverse Hankel transform when t-->infty
    y=besselj(1.0,kap*R0).*besselj(0.0,kap*rk)./(rg + D*kap.^4);
    end

end
