function M = calc_M_E1(p, phys, sys, err, N)
% Nonrelativistic, one-sided axial-gradient translation for the revision.
assert(~isfield(phys,'is_relativistic') || ~phys.is_relativistic, ...
    'This revision override is restricted to the nonrelativistic map.');
M=repmat(eye(3),1,1,N);
E_average=zeros(1,N);
if isfield(err,'func_dEz') && isfield(err,'z1') && isfield(err,'z2')
    z_nodes=linspace(err.z1,err.z2,100);
    for i=1:numel(z_nodes)
        E_average=E_average+err.func_dEz(phys.X,phys.Y,z_nodes(i))/100;
    end
end
e=1.602e-19; m=9.109e-31;
Bz=phys.B_vec(3,:);
delta_L=2*pi^2*sys.N_turn^2*m./(e*Bz.^2).*E_average;
delta=-delta_L.*phys.tan_gamma;
M(1,3,:)=delta.*cos(p.phi);
M(2,3,:)=delta.*sin(p.phi);
end
