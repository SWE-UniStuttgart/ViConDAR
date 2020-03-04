%% WriteBladedWindFile
% Function: Writes  Bladed wind field from windfield struct.
% Based on "Bladed wind file format" Version 3.67
%
%% Usage:
%
% WriteBladedWindFile(windfield,Parameter)
%
%% Input:
% windfield - struct
% parameter - struct
%
%% Output:
%
%
%% Modified:
% * Tim Hagemann on 11-Oct-2016
% - corrected Dummy to windfield.Fmax
% - added .$02 support for fourier file
% - changed substraction of turbulence and shear to revert
% ReadBladedWindFile correctly
% - StaticWindfield calculation with 3D Matrix
%
%% ToDo:
%
%
%
%% Created:
% David Schlipf on 13-Oct-2014
%
% (c) Universitaet Stuttgart
%
%% Code
function WriteBladedWindFile(windfield,Parameter)

FileName                        = Parameter.Wind.FileName;

FileFormat                      = -99;  % current format
SpectraTyp                      = 7;    % Kaimal
nHeader                         = 92;
nTurbulenceComponent            = 3;

if ~isfield(windfield,'Fmax')
    windfield.Fmax              = windfield.URef/(2*windfield.grid.dx);
end

fid                             = fopen(FileName,'w');
fwrite( fid, FileFormat,            'int16' );
fwrite( fid, SpectraTyp,            'int16' );
fwrite( fid, nHeader,               'int32' );
fwrite( fid, nTurbulenceComponent,  'int32' );

fwrite( fid, windfield.grid.dz,     'real*4');
fwrite( fid, windfield.grid.dy,     'real*4');
fwrite( fid, windfield.grid.dx,     'real*4');
fwrite( fid, windfield.nFFT/2,      'int32' );
fwrite( fid, windfield.URef,        'real*4');
fwrite( fid, windfield.L(1,3),      'real*4');
fwrite( fid, windfield.L(1,2),      'real*4');
fwrite( fid, windfield.L(1,1),      'real*4');
fwrite( fid, windfield.Fmax,        'real*4');
fwrite( fid, windfield.RandSeed,    'int32' );
fwrite( fid, windfield.grid.nz,     'int32' );
fwrite( fid, windfield.grid.ny,     'int32' );

fwrite( fid, windfield.L(2,3),      'real*4');
fwrite( fid, windfield.L(2,2),      'real*4');
fwrite( fid, windfield.L(2,1),      'real*4');
fwrite( fid, windfield.L(3,3),      'real*4');
fwrite( fid, windfield.L(3,2),      'real*4');
fwrite( fid, windfield.L(3,1),      'real*4');

fwrite( fid, windfield.c_c,         'real*4');
fwrite( fid, windfield.L_c,         'real*4');



% revert ReadBladedWindFile's calculation
HubHeight           = Parameter.Turbine.HubHeight;
WindShearExponent   = Parameter.Wind.WindShearExponent;
[T3,Y3,Z3]          = meshgrid(windfield.grid.t,windfield.grid.y,windfield.grid.z);
StaticWindfield     = ((Z3+HubHeight)/HubHeight).^WindShearExponent*windfield.URef;

U_n                 = round((windfield.u-StaticWindfield)/(Parameter.Wind.I(1)*windfield.URef)*1000);
V_n                 = round((windfield.v-0              )/(Parameter.Wind.I(2)*windfield.URef)*1000);
W_n                 = round((windfield.w-0              )/(Parameter.Wind.I(3)*windfield.URef)*1000);



Data            = NaN(windfield.grid.nz*windfield.grid.nz*windfield.nFFT*3,1);
Data(1:3:end)   = reshape(permute(U_n,[1 3 2]),1,[]);
Data(2:3:end)   = reshape(permute(V_n,[1 3 2]),1,[]);
Data(3:3:end)   = reshape(permute(W_n,[1 3 2]),1,[]);


fwrite( fid, Data , 'int16');

fclose(fid);


if isfield(Parameter.Wind,'FourierFileFlag')
    if Parameter.Wind.FourierFileFlag
        % write .$02 file
        fid = fopen([FileName(1:end-4) '.$02'],'w');
        gatheredData(1,:,:,:,:) = windfield.FFT.Amplitude;
        gatheredData(2,:,:,:,:) = windfield.FFT.Phase;
        FFTData = reshape(gatheredData,1,[]);
        fwrite(fid,FFTData,'real*4');
        fclose(fid);
        
        % write comment file
        fid  = fopen([FileName(1:end-4) '.%02'],'w');
        fprintf(fid,'FILE   %s\n',[FileName(1:end-4) '.$02']);
        fprintf(fid,'ACCESS    D\n');
        fprintf(fid,'FORM      U\n');
        fprintf(fid,'RECL                 4\n');
        fprintf(fid,'FORMAT    R4\n');
        fprintf(fid,'HEADREC              0\n');
        fprintf(fid,'CONTENT   WIND\n');
        fprintf(fid,'CONFIG    STATIONARY\n');
        fprintf(fid,'NDIMENS              5\n');
        fprintf(fid,'DIMENS               2        %i          %i          %i           %i\n',windfield.nFFT/2,windfield.grid.ny,windfield.grid.nz,nTurbulenceComponent);
        fprintf(fid,'GENLAB    ''Normalised spectrum and phase''\n');
        fprintf(fid,'VARIAB    Spectrum ''Phase angle''\n');
        fprintf(fid,'VARUNIT   T A\n');
        fprintf(fid,'AXISLAB   ''Frequency''\n');
        fprintf(fid,'AXIUNIT   1/T\n');
        fprintf(fid,'AXIMETH              2\n');
        fprintf(fid,'MIN             %.6f\n',windfield.FFT.f(1));
        fprintf(fid,'STEP            %.6f\n',windfield.FFT.f(1));
        fprintf(fid,'AXISLAB   ''Lateral position''\n');
        fprintf(fid,'AXIUNIT   L\n');
        fprintf(fid,'AXIMETH              2\n');
        fprintf(fid,'MIN          %.6f\n',min(windfield.grid.y));
        fprintf(fid,'STEP           %.6f\n',windfield.grid.dy);
        fprintf(fid,'AXISLAB   ''Vertical position''\n');
        fprintf(fid,'AXIUNIT   L\n');
        fprintf(fid,'AXIMETH              2\n');
        fprintf(fid,'MIN          %.6f\n',min(windfield.grid.z));
        fprintf(fid,'STEP           %.6f\n',windfield.grid.dz);
        fprintf(fid,'AXISLAB   ''Component''\n');
        fprintf(fid,'AXIUNIT   N\n');
        fprintf(fid,'AXIMETH              1\n');
        fprintf(fid,'AXITICK   Longitudinal Lateral      Vertical\n');
        fprintf(fid,'NVARS                3\n');
        fprintf(fid,'XLU           %.6f\n',windfield.L(1,1));
        fprintf(fid,'XLV           %.6f\n',windfield.L(2,1));
        fprintf(fid,'XLW            %.6f\n',windfield.L(3,1));
        fclose(fid);
        
        % create empty pj file
        fclose( fopen([FileName(1:end-4) '.$PJ'],'w'));
    end
end

end
