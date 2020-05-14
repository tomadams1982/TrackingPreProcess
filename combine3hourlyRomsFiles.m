function [data] = combine3hourlyRomsFiles(romsDir,date,varargin)
%COMBINEROMSFILES Combine multiple hourly ROMS files into a single
% struct for a day, which contains values aggregated over 3 depth ranges,
% and uses temporal interpolation to fill in gaps between 3 hourly source
% files.
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
    
    prefix='NEATL_SAMS_';
    
    %date='20190601';
    %date=datetime(date,'InputFormat','yyyyMMdd');
    d1=datestr(date,'yyyymmdd');
    
    date2=date+1;
    d2=datestr(date2,'yyyymmdd');
    %d2='20190602';
    
    times={'00','03','06','09','12','15','18','21'};
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
        
%         data.u(:,:,3,str2num(times{tt})+1)=mean(u(:,:,1:10,1),3); % Bed
%         data.u(:,:,2,str2num(times{tt})+1)=mean(u(:,:,11:30,1),3); % Mid
%         data.u(:,:,1,str2num(times{tt})+1)=mean(u(:,:,31:40,1),3); % Surface

        v=ncread(file,'v');
        data.v(:,:,4,times2(tt)+1)=mean(v(:,:,1:10,1),3); % Bed
        data.v(:,:,3,times2(tt)+1)=mean(v(:,:,11:30,1),3); % Mid
        data.v(:,:,2,times2(tt)+1)=mean(v(:,:,31:36,1),3); % Near-surface
        data.v(:,:,1,times2(tt)+1)=mean(v(:,:,37:40,1),3); % Surface
        
%         data.v(:,:,3,str2num(times{tt})+1)=mean(v(:,:,1:10,1),3); % Bed
%         data.v(:,:,2,str2num(times{tt})+1)=mean(v(:,:,11:30,1),3); % Mid
%         data.v(:,:,1,str2num(times{tt})+1)=mean(v(:,:,31:40,1),3); % Surface

        if (velocityOnly==0)
            s=ncread(file,'salt');
            data.s(:,:,4,times2(tt)+1)=mean(s(:,:,1:10,1),3); % Bed
            data.s(:,:,3,times2(tt)+1)=mean(s(:,:,11:30,1),3); % Mid
            data.s(:,:,2,times2(tt)+1)=mean(s(:,:,31:36,1),3); % Near-surface
            data.s(:,:,1,times2(tt)+1)=mean(s(:,:,37:40,1),3); % Surface
            
%             data.s(:,:,3,str2num(times{tt})+1)=mean(s(:,:,1:10,1),3); % Bed
%             data.s(:,:,2,str2num(times{tt})+1)=mean(s(:,:,11:30,1),3); % Mid
%             data.s(:,:,1,str2num(times{tt})+1)=mean(s(:,:,31:40,1),3); % Surface 

            t=ncread(file,'temp');
            data.t(:,:,4,times2(tt)+1)=mean(t(:,:,1:10,1),3); % Bed
            data.t(:,:,3,times2(tt)+1)=mean(t(:,:,11:30,1),3); % Mid
            data.t(:,:,2,times2(tt)+1)=mean(t(:,:,31:36,1),3); % Near-surface
            data.t(:,:,1,times2(tt)+1)=mean(t(:,:,37:40,1),3); % Surface
        
%             data.t(:,:,3,str2num(times{tt})+1)=mean(t(:,:,1:10,1),3); % Bed
%             data.t(:,:,2,str2num(times{tt})+1)=mean(t(:,:,11:30,1),3); % Mid
%             data.t(:,:,1,str2num(times{tt})+1)=mean(t(:,:,31:40,1),3); % Surface 
        end
    end

    % Do linear temporal interpolation
    for tt=0:6
        data.u(:,:,:,2+(tt*3))=(2/3)*data.u(:,:,:,1+(tt*3)) + (1/3)*data.u(:,:,:,4+(tt*3));
        data.u(:,:,:,3+(tt*3))=(1/3)*data.u(:,:,:,1+(tt*3)) + (2/3)*data.u(:,:,:,4+(tt*3));

        data.v(:,:,:,2+(tt*3))=(2/3)*data.v(:,:,:,1+(tt*3)) + (1/3)*data.v(:,:,:,4+(tt*3));
        data.v(:,:,:,3+(tt*3))=(1/3)*data.v(:,:,:,1+(tt*3)) + (2/3)*data.v(:,:,:,4+(tt*3));

        if (velocityOnly==0)
            data.s(:,:,:,2+(tt*3))=(2/3)*data.s(:,:,:,1+(tt*3)) + (1/3)*data.s(:,:,:,4+(tt*3));
            data.s(:,:,:,3+(tt*3))=(1/3)*data.s(:,:,:,1+(tt*3)) + (2/3)*data.s(:,:,:,4+(tt*3));

            data.t(:,:,:,2+(tt*3))=(2/3)*data.t(:,:,:,1+(tt*3)) + (1/3)*data.t(:,:,:,4+(tt*3));
            data.t(:,:,:,3+(tt*3))=(1/3)*data.t(:,:,:,1+(tt*3)) + (2/3)*data.t(:,:,:,4+(tt*3));
        end
    end

    % Stuff using the first value from the next day
    u2av=zeros(size(data.lon_u,1),size(data.lon_u,2),3,24);
    v2av=zeros(size(data.lon_v,1),size(data.lon_v,2),3,24);
    if (velocityOnly==0)
        s2av=zeros(size(data.lon_rho,1),size(data.lon_rho,2),3,24);
        t2av=zeros(size(data.lon_rho,1),size(data.lon_rho,2),3,24);
    end
    
    file2=[romsDir prefix d2 '00.nc'];
    u2=ncread(file2,'u');
    u2av(:,:,4,times2(tt)+1)=mean(u2(:,:,1:10,1),3); % Bed
    u2av(:,:,3,times2(tt)+1)=mean(u2(:,:,11:30,1),3); % Mid
    u2av(:,:,2,times2(tt)+1)=mean(u2(:,:,31:36,1),3); % Near-surface
    u2av(:,:,1,times2(tt)+1)=mean(u2(:,:,37:40,1),3); % Surface    
%     u2av(:,:,3,str2num(times{tt})+1)=mean(u2(:,:,1:10,1),3); % Bed
%     u2av(:,:,2,str2num(times{tt})+1)=mean(u2(:,:,11:30,1),3); % Mid
%     u2av(:,:,1,str2num(times{tt})+1)=mean(u2(:,:,31:40,1),3); % Surface
    v2=ncread(file2,'v');
    v2av(:,:,4,times2(tt)+1)=mean(v2(:,:,1:10,1),3); % Bed
    v2av(:,:,3,times2(tt)+1)=mean(v2(:,:,11:30,1),3); % Mid
    v2av(:,:,2,times2(tt)+1)=mean(v2(:,:,31:36,1),3); % Near-surface
    v2av(:,:,1,times2(tt)+1)=mean(v2(:,:,37:40,1),3); % Surface  
%     v2av(:,:,3,str2num(times{tt})+1)=mean(v2(:,:,1:10,1),3); % Bed
%     v2av(:,:,2,str2num(times{tt})+1)=mean(v2(:,:,11:30,1),3); % Mid
%     v2av(:,:,1,str2num(times{tt})+1)=mean(v2(:,:,31:40,1),3); % Surface
    if (velocityOnly==0)
        s2=ncread(file2,'salt');
        s2av(:,:,4,times2(tt)+1)=mean(s2(:,:,1:10,1),3); % Bed
        s2av(:,:,3,times2(tt)+1)=mean(s2(:,:,11:30,1),3); % Mid
        s2av(:,:,2,times2(tt)+1)=mean(s2(:,:,31:36,1),3); % Near-surface
        s2av(:,:,1,times2(tt)+1)=mean(s2(:,:,37:40,1),3); % Surface  
%         s2av(:,:,3,str2num(times{tt})+1)=mean(s2(:,:,1:10,1),3); % Bed
%         s2av(:,:,2,str2num(times{tt})+1)=mean(s2(:,:,11:30,1),3); % Mid
%         s2av(:,:,1,str2num(times{tt})+1)=mean(s2(:,:,31:40,1),3); % Surface
        t2=ncread(file2,'temp');
        t2av(:,:,4,times2(tt)+1)=mean(t2(:,:,1:10,1),3); % Bed
        t2av(:,:,3,times2(tt)+1)=mean(t2(:,:,11:30,1),3); % Mid
        t2av(:,:,2,times2(tt)+1)=mean(t2(:,:,31:36,1),3); % Near-surface
        t2av(:,:,1,times2(tt)+1)=mean(t2(:,:,37:40,1),3); % Surface  
%         t2av(:,:,3,str2num(times{tt})+1)=mean(t2(:,:,1:10,1),3); % Bed
%         t2av(:,:,2,str2num(times{tt})+1)=mean(t2(:,:,11:30,1),3); % Mid
%         t2av(:,:,1,str2num(times{tt})+1)=mean(t2(:,:,31:40,1),3); % Surface
    end
    
    data.u(:,:,:,23)=(2/3)*data.u(:,:,:,22) + (1/3)*u2av(:,:,:,1);
    data.u(:,:,:,24)=(1/3)*data.u(:,:,:,22) + (2/3)*u2av(:,:,:,1);

    data.v(:,:,:,23)=(2/3)*data.v(:,:,:,22) + (1/3)*v2av(:,:,:,1);
    data.v(:,:,:,24)=(1/3)*data.v(:,:,:,22) + (2/3)*v2av(:,:,:,1);

    if (velocityOnly==0)
        data.s(:,:,:,23)=(2/3)*data.s(:,:,:,22) + (1/3)*s2av(:,:,:,1);
        data.s(:,:,:,24)=(1/3)*data.s(:,:,:,22) + (2/3)*s2av(:,:,:,1);

        data.t(:,:,:,23)=(2/3)*data.t(:,:,:,22) + (1/3)*t2av(:,:,:,1);
        data.t(:,:,:,24)=(1/3)*data.t(:,:,:,22) + (2/3)*t2av(:,:,:,1);
    end
    % Save file
    %save([romsDir prefix date '.mat'],'data')
    
end

