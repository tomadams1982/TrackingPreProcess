% processROMSFiles
basedir = 'C:\Users\sa01ta\Documents\COMPASS\scripts\'
cd(basedir)
% MeshCrop = load([basedir 'Mesh_NEATL_crop.mat']);
% MeshCrop = MeshCrop.MeshCrop;

% Are the files to be processed the hourly ones, or the three-hourly ones?
onehour=0;
if (onehour == 1)
    load([basedir 'Mesh_NEATL_crop2_4layers.mat']);
    Mesh=MeshCrop2;
    romsDir='D:\hydroOut\NEA_ROMS\1hr\2019\';
    outDir='D:\hydroOut\NEA_ROMS\24hr4layerTri\2019\';
    prefix='NEATL_';
else
    load('Mesh_NEATL_SAMS.mat')
    romsDir='D:\hydroOut\NEA_ROMS\SAMS_CROP\3hr\2019\';
    outDir='D:\hydroOut\NEA_ROMS\SAMS_CROP\24hr4layerTri\2019\';
    prefix='NEATL_';
end

startDateString='20190501';
nDays=31;
startDate=datetime(startDateString,'InputFormat','yyyyMMdd');
dates=startDate+(0:(nDays-1));

for day=1:nDays
    disp(['**** Processing ROMS files for date ' datestr(dates(day),'yyyymmdd') ' ****']);
    
    clear data hydroDat
    
    % Combining 1 hourly, larger domain files into smaller structs
    % Takes about 110 s per day
    tic
    
    if (onehour==1)
        disp('- Combining 1-hourly files to single data structure');
        data=combine1hourlyRomsFiles(romsDir,dates(day),'velocityOnly',1);
    else
        disp('- Combining 3-hourly files to single data structure');
        data=combine3hourlyRomsFiles(romsDir,dates(day),'velocityOnly',0);
    end
    %save([outDir prefix datestr(dates(day),'yyyymmdd') '.mat'],'data');
    toc
    
    % Reinterpolate on to new mesh - this will incorporate the cropping required
    % Takes 500 seconds per day
    %size(data.u)
    tic
    disp('- Interpolating values to specified triangular mesh');
    hydroDat = interpolateHydro(Mesh,data);
    toc
    
    % Create .nc file
    % Takes just a few seconds
    tic
    disp('- Writing interpolated structure to NetCDF format');
    delete([outDir prefix datestr(dates(day),'yyyymmdd') '.nc']);
    writeHydroNC([outDir prefix datestr(dates(day),'yyyymmdd') '.nc'],Mesh,hydroDat);
    toc
    
end


