%% Header
%
% Function to create the windfield global format
% Forked from NREL OpenFAST repository and modified by SWE (2019)

function [windfield] = velocity2windfield(velocity,dz,dy,dt,SummVars)

% transform windfield structure (SWE)
%        to velocity 4-D vector (NREL)
%
% inputs:   - windfield: .grid.nt
%                             .ny
%                             .nz
%                        .u
%                        .v
%                        .w
% outputs:  - velocity: 4-D vector: (nt x (u,v,w) x ny x nz)
%           - dz, dy, dt    - scalars: distance between two points in the vertical [m]/
%                                      horizontal [m]/time [s] dimension
%
%
% for use with
% function writeBLgrid(FileName, velocity, dz, dy, dt, zOffset, z0, SummVars)
% - SumVars: 6 variables from the summary file {zHub, Clockwise, UBAR, TI_u, TI_v, TI_w}
%                                                 90          1    12    12   9.6     6
%
% missing:  - zOffset
%           - z0
%           - SummVars
%
% (run duration for nt=2561, ny=nz=23: 1.6 seconds)
%
[nt , ~, ny, nz]    = size(velocity);
windfield.grid.nt   = nt;
windfield.grid.ny   = ny;
windfield.grid.nz   = nz;

windfield.ny        = ny;
windfield.nz        = nz;

windfield.grid.dt   = round(dt*100000)/100000; % issue with constrained time steps off
windfield.dt        = round(dt*100000)/100000; % issue with constrained time steps off
windfield.grid.dy   = dy;
windfield.grid.dz   = dz;

windfield.grid.t    = [0:windfield.dt:windfield.dt*(windfield.grid.nt-1)]';

windfield.grid.z    = -windfield.grid.dz*(windfield.grid.nz-1)/2:windfield.grid.dz:windfield.grid.dz*(windfield.grid.nz-1)/2;
windfield.grid.y    = -windfield.grid.dy*(windfield.grid.ny-1)/2:windfield.grid.dy:windfield.grid.dy*(windfield.grid.ny-1)/2;

[windfield.grid.Y, windfield.grid.Z] = meshgrid(windfield.grid.y,windfield.grid.z);

windfield.grid.dy   = dy;
windfield.grid.dz   = dz;

windfield.u         = zeros(ny,nt,nz);
windfield.v         = zeros(ny,nt,nz);
windfield.w         = zeros(ny,nt,nz);
for it=1:nt
    windfield.u(:,it,:) = squeeze(velocity(it,1,:,:));
    windfield.v(:,it,:) = squeeze(velocity(it,2,:,:));
    windfield.w(:,it,:) = squeeze(velocity(it,3,:,:));   
end

windfield.URef      = SummVars(3);    
windfield.T_offset  = windfield.grid.dy*(windfield.grid.ny-1)/windfield.URef/2+windfield.grid.dt; %GridWidth/URef/2 % +dt(DS)
