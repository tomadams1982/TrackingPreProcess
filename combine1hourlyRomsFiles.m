function [data] = combine1hourlyRomsFiles(romsDir,date,varargin)
%COMBINE1HOURLYROMSFILES Combine multiple hourly ROMS files into a single
% struct for a day, which contains values aggregated over 3 depth ranges.
%
% Inputs:   dataDirectory   - directory containing the ROMS files
%           date            - date to generate the file for (datetime

    % Files copied 
    % From: (samhanach) /home/sa01da/data/oss
    % To: W:\sa01ta\hydroFilesTemp\OSS

    %file='F:\hydroOut\NEA_ROMS\2017_OLD\NEATL_SAMS_2017123000.nc';
    %romsDir='W:\sa01ta\hydroFilesTemp\OSS\2019\';

    velocityOnly = 0;
    
    for i = 1:2:length(varargin) % only bother with odd arguments, i.e. the labels
        switch varargin{i}
            case 'velocityOnly'
                velocityOnly = varargin{i+1};
        end
    end
    
    prefix='NEATL_';
    
    %date='20190601';
    %date=datetime(date,'InputFormat','yyyyMMdd');
    d1=datestr(date,'yyyymmdd');
    
    %date2=date+1;
    %d2=datestr(date2,'yyyymmdd');
    %d2='20190602';
    
    %times={'00','03','06','09','12','15','18','21'};
    times=strsplit(sprintf('%02d ',0:23),' ');
    times=times(1,1:24);
    times2=cellfun(@str2num,times);
    
    assert(length(dir([romsDir prefix d1 '*.nc']))==length(times),'Not all files for requested day were present');
    
    file=[romsDir prefix d1 times{1} '.nc'];
    %ncdisp(file)

    % Grid data
    data.lon_u=ncread(file,'lon_u');
    data.lat_u=ncread(file,'lat_u');
    data.lon_v=ncread(file,'lon_v');
    data.lat_v=ncread(file,'lat_v');
    data.lon_rho=ncread(file,'lon_rho');
    data.lat_rho=ncread(file,'lat_rho');
    data.s_rho=ncread(file,'s_rho');
    data.h=ncread(file,'h');

    data.u=zeros(size(data.lon_u,1),size(data.lon_u,2),3,24);
    data.v=zeros(size(data.lon_v,1),size(data.lon_v,2),3,24);
    if (velocityOnly==0)
        data.s=zeros(size(data.lon_rho,1),size(data.lon_rho,2),3,24);
        data.t=zeros(size(data.lon_rho,1),size(data.lon_rho,2),3,24);
    end
    
    % Get depth averaged values
    for tt=1:length(times)
        file=[romsDir prefix d1 times{tt} '.nc'];
        disp(file);

        u=ncread(file,'u');
        data.u(:,:,4,times2(tt)+1)=mean(u(:,:,1:10,1),3); % Bed
        data.u(:,:,3,times2(tt)+1)=mean(u(:,:,11:30,1),3); % Mid
        data.u(:,:,2,times2(tt)+1)=mean(u(:,:,31:36,1),3); % Near-surface
        data.u(:,:,1,times2(tt)+1)=mean(u(:,:,37:40,1),3); % Surface

        v=ncread(file,'v');
        data.v(:,:,4,times2(tt)+1)=mean(v(:,:,1:10,1),3); % Bed
        data.v(:,:,3,times2(tt)+1)=mean(v(:,:,11:30,1),3); % Mid
        data.v(:,:,2,times2(tt)+1)=mean(v(:,:,31:36,1),3); % Near-surface
        data.v(:,:,1,times2(tt)+1)=mean(v(:,:,37:40,1),3); % Surface

        if (velocityOnly==0)
            try
                s=ncread(file,'salt');
                data.s(:,:,4,times2(tt)+1)=mean(s(:,:,1:10,1),3); % Bed
                data.s(:,:,3,times2(tt)+1)=mean(s(:,:,11:30,1),3); % Mid
                data.s(:,:,2,times2(tt)+1)=mean(s(:,:,31:36,1),3); % Near-surface 
                data.s(:,:,1,times2(tt)+1)=mean(s(:,:,37:40,1),3); % Surface 
            catch
                warning('Salinity information not present in file');
            end

            try
                t=ncread(file,'temp');
                data.t(:,:,4,times2(tt)+1)=mean(t(:,:,1:10,1),3); % Bed
                data.t(:,:,3,times2(tt)+1)=mean(t(:,:,11:30,1),3); % Mid
                data.t(:,:,2,times2(tt)+1)=mean(t(:,:,31:36,1),3); % Near-surface 
                data.t(:,:,1,times2(tt)+1)=mean(t(:,:,37:40,1),3); % Surface 
            catch
                warning('Temperature information not present in file');
            end
        end
    end

    % Save file
    %save([romsDir prefix date '.mat'],'data')
    
end

