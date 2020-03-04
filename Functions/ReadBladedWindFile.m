%% ReadBladedWindFile
% Function: Reads in Bladed wind field to windfield struct.
% Based on "Bladed wind file format" Version 3.67
%
%% Usage:
%
% windfield  = ReadBladedWindFile(Parameter)
%
%% Input:
%
% Parameter.Wind.WindShearExponent  = 0.16;     %[-]
% Parameter.Turbine.HubHeight       = 119.03;   %[m]
% Parameter.Wind.I                  = [1 0.8 0.5]*17.6/100;    % [-]
% Parameter.Wind.FileName           = 'D:\SheffieldMeeting\Winds\16_224_224_660.wnd';
%
%% Output:
% windfield - struct
%
%
%% Modified:
% * Tim Hagemann on 11-Oct-2016
% - corrected Dummy to windfield.Fmax
% - added support for .$02 file: FourierFileFlag
% - StaticWindfield calculation with 3D Matrix
%
%% ToDo:
%
%
%
%% Created:
% David Schlipf on 09-May-2014
%
% (c) Universitaet Stuttgart
%
%% Code
function windfield  = ReadBladedWindFile(Parameter)

WindShearExponent           = Parameter.Wind.WindShearExponent;
HubHeight                   = Parameter.Turbine.HubHeight;
I                           = Parameter.Wind.I;
FileName                    = Parameter.Wind.FileName;

% Check Ending
Ending              = FileName(end-3:end);
if ~strcmpi( Ending, '.wnd' )
    error('No *.wnd file!')
end

fid                 = fopen(FileName);

FileFormat                      = fread( fid, 1, 'int16' );

if FileFormat ~= -99
    error('Wrong Format!')
end

SpectraTyp                      = fread( fid, 1, 'int16' );

if SpectraTyp >= 2
    nHeader                     = fread( fid, 1, 'int32' );
    nTurbulenceComponent        = fread( fid, 1, 'int32' );
end

if SpectraTyp == 4
    nTurbulenceComponentIvK     = fread( fid, 1, 'int16' );
    ImprovedvanKarmanModel      = fread( fid, 6, 'real*4');
end

windfield.grid.dz               = fread( fid, 1, 'real*4');
windfield.grid.dy               = fread( fid, 1, 'real*4');
windfield.grid.dx               = fread( fid, 1, 'real*4');

windfield.nFFT                  = fread( fid, 1, 'int32' )*2;
windfield.URef                  = fread( fid, 1, 'real*4');
windfield.L(1,3)                = fread( fid, 1, 'real*4');
windfield.L(1,2)                = fread( fid, 1, 'real*4');
windfield.L(1,1)                = fread( fid, 1, 'real*4');

windfield.Fmax                  = fread( fid, 1, 'real*4');
windfield.RandSeed              = fread( fid, 1, 'int32' );
windfield.grid.nz               = fread( fid, 1, 'int32' );
windfield.grid.ny               = fread( fid, 1, 'int32' );

if nTurbulenceComponent == 3
    windfield.L(2,3)         	= fread( fid, 1, 'real*4');
    windfield.L(2,2)           	= fread( fid, 1, 'real*4');
    windfield.L(2,1)           	= fread( fid, 1, 'real*4');
    windfield.L(3,3)         	= fread( fid, 1, 'real*4');
    windfield.L(3,2)           	= fread( fid, 1, 'real*4');
    windfield.L(3,1)            = fread( fid, 1, 'real*4');
end

if SpectraTyp == 7
    windfield.c_c               = fread( fid, 1, 'real*4');
    windfield.L_c               = fread( fid, 1, 'real*4');
end

if SpectraTyp == 8
    MannModel{1}                = fread( fid, 6, 'real*4');
    MannModel{2}                = fread( fid, 3, 'int32' );
    MannModel{3}                = fread( fid, 2, 'real*4');
    MannModel{4}                = fread( fid, 3, 'int32' );
    MannModel{5}                = fread( fid, 2, 'real*4');
end

nDataPoints         = windfield.grid.nz*windfield.grid.nz*windfield.nFFT*3;
Data                = fread( fid, nDataPoints, 'int16')/1000;

windfield.grid.t    = [0:windfield.nFFT-1]*windfield.grid.dx/windfield.URef;
windfield.grid.y    = [-(windfield.grid.ny-1)/2:(windfield.grid.ny-1)/2]*windfield.grid.dy;
windfield.grid.z    = [-(windfield.grid.nz-1)/2:(windfield.grid.nz-1)/2]*windfield.grid.dz;
[windfield.grid.Y, windfield.grid.Z]    = meshgrid(windfield.grid.y,windfield.grid.z);

[T3,Y3,Z3]          = meshgrid(windfield.grid.t,windfield.grid.y,windfield.grid.z);

StaticWindfield     = ((Z3+HubHeight)/HubHeight).^WindShearExponent*windfield.URef;
windfield.u         = permute(reshape(Data(1:3:end)*windfield.URef*I(1),windfield.grid.ny,windfield.grid.nz,windfield.nFFT),            [1 3 2])+StaticWindfield;
windfield.v         = permute(reshape(Data(2:3:end)*windfield.URef*I(2),windfield.grid.ny,windfield.grid.nz,windfield.nFFT),            [1 3 2]);
windfield.w         = permute(reshape(Data(3:3:end)*windfield.URef*I(3),windfield.grid.ny,windfield.grid.nz,windfield.nFFT),            [1 3 2]);

fclose(fid);



if isfield(Parameter.Wind,'FourierFileFlag')
    if Parameter.Wind.FourierFileFlag
        FTFileName              = [Parameter.Wind.FileName(1:end-4),'.$02'];
        fid                     = fopen(FTFileName);
        FFTData                 = reshape(fread(fid,'real*4'),[2,windfield.nFFT/2,windfield.grid.ny,windfield.grid.nz,nTurbulenceComponent]);
        windfield.FFT.Amplitude = squeeze(FFTData(1,:,:,:,:));
        windfield.FFT.Phase     = squeeze(FFTData(2,:,:,:,:));
        windfield.FFT.f         = windfield.Fmax*2/windfield.nFFT*(1:windfield.nFFT/2);
        fclose(fid);
    end
end
%% Code to verify
% DataASCII   	= importdata([FileName(1:end-4),'_vH.txt']);
% DataASCII   	= importdata([FileName(1:end-4),'AtBottomRight.txt']);
% u_ASCII     = DataASCII.data(:,2);
% t_ASCII     = DataASCII.data(:,1);
% u_Binary    = (windfield.u(34,:,1));
% t_Binary    = [0:windfield.nFFT*2-1]*windfield.grid.dx/windfield.URef;
%
% figure
% hold on;box on;grid on
% plot(t_Binary,u_Binary,'.-r')
% plot(50+t_ASCII+(178/2)/windfield.URef,u_ASCII,'o-b')
% xlim([60 90])
