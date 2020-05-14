function [Mesh] = createRomsTriangularMesh(romsDir, date, prefix, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    narginchk(3,7);
    
    xl=[];
    yl=[];
    for i = 1:2:length(varargin) % only bother with odd arguments, i.e. the labels
        switch varargin{i}
            case 'xl'
                xl = varargin{i+1};
            case 'yl'
                yl = varargin{i+1};
        end
    end

    %%
    % Create a triangular mesh based on ROMS rho points
    % In Dima's FVCOM grid, triangles are based on a set of nodes defined
    % anti-clockwise. Match this here
    %prefix='NEATL_';
    file=[romsDir prefix datestr(date,'yyyymmdd') '00.nc'];
    
    %file=file2;
    
    disp('Reading rho points')
    
    lon_rho=ncread(file,'lon_rho');
    lat_rho=ncread(file,'lat_rho');
    depth_rho=ncread(file,'h');
    
%     size(lon_rho)
%     size(lat_rho)
%     size(depth_rho)
    
    mask_rho=ncread(file,'mask_rho');
    
    if ~isempty(xl)
        lon_rho=lon_rho(xl(1):xl(2),:);
        lat_rho=lat_rho(xl(1):xl(2),:);
        mask_rho=mask_rho(xl(1):xl(2),:);
        depth_rho=depth_rho(xl(1):xl(2),:);
    end
    if ~isempty(yl)
        lon_rho=lon_rho(:,yl(1):yl(2));
        lat_rho=lat_rho(:,yl(1):yl(2));
        mask_rho=mask_rho(:,yl(1):yl(2));
        depth_rho=depth_rho(:,yl(1):yl(2));
    end
    
    % Assign rho points a unique index = i+(j-1)*n1
    n1=size(lon_rho,1);
    n2=size(lon_rho,2);
    
    tic
    disp('Creating trinodes (possibly not necessary as will redo in a bit)');
    trinodesLowerLeft=zeros((n1-1)*(n2-1),3);
    trinodesUpperRight=zeros((n1-1)*(n2-1),3);
    for j=0:(n2-2)
        for i=1:(n1-1)
            trinodesLowerLeft(j*n1+i,:)=[j*n1+i j*n1+i+1 (j+1)*n1+i+1]; % "lower" triangles
            trinodesUpperRight(j*n1+i,:)=[j*n1+i (j+1)*n1+i+1 (j+1)*n1+i]; % "upper" triangles

            %trinodes2((j-1)*n1+i,:)=[i j];
        end
    end
    trinodesLowerLeft(n1:n1:end,:)=[];
    trinodesUpperRight(n1:n1:end,:)=[];

    trinodes=zeros(length(trinodesLowerLeft)+length(trinodesUpperRight),3);
    % Interleaved - easier to derive indices for neighbouring elements
    trinodes(1:2:(length(trinodesLowerLeft)*2),:)=trinodesLowerLeft;
    trinodes(2:2:(length(trinodesUpperRight)*2),:)=trinodesUpperRight;
    % Sequential
    %trinodes=[trinodesLowerLeft;trinodesUpperRight];
    toc
    
    tic
    disp('Creating corner nodes');
    % reshape rho points to match index
    lon_rho2 = reshape(lon_rho,[],1);
    lat_rho2 = reshape(lat_rho,[],1);

    % Coordinates of all corner nodes
    nodexy=[lon_rho2,lat_rho2];

    % Coordinates of centroids
    uvnode=zeros(length(trinodes),2);
    for i=1:length(uvnode)
        uvnode(i,1)=mean(nodexy(trinodes(i,:),1));
        uvnode(i,2)=mean(nodexy(trinodes(i,:),2));
    end
    toc

    tic
    disp('Creating neighbouring element list');
    % Neighbouring elements
    % - labelled anticlockwise starting from the diagonal side
    nbe=zeros(length(trinodes),3);
    % Lower left triangles
    for i=1:2:length(nbe)
        nbe(i,1)=i+1; % top right (diagonal) edge
        if i>2*(n1-1)
            nbe(i,2)=i-2*(n1-1)+1; % left (vertical) edge
        end
        if mod(i,2*(n1-1))~=(2*(n1-1)-1)
            nbe(i,3)=i+3; % bottom (horizontal) edge
        end
    end
    % Upper right triangles
    for i=2:2:length(nbe)
        nbe(i,1)=i-1; % bottom left (diagonal) edge
        if i<(n1-1)*2*(n2-2)
            nbe(i,2)=i+2*(n1-1)-1; % right (vertical) edge
        end
        if mod(i,2*(n1-1))~=2
            nbe(i,3)=i-3; % top (horizontal) edge
        end
    end
    toc

    % Let's have a look
    % clf
    % scatter(uvnode(1,1),uvnode(1,2),'b')
    % hold on
    % plot(nodexy(trinodes(1,[1:end 1]),1),nodexy(trinodes(1,[1:end 1]),2),'b')
    % scatter(uvnode(nbe(1,[1 3]),1),uvnode(nbe(1,[1 3]),2),'k')
    % plot(nodexy(trinodes(2,[1:end 1]),1),nodexy(trinodes(2,[1:end 1]),2),'k')
    % plot(nodexy(trinodes(4,[1:end 1]),1),nodexy(trinodes(4,[1:end 1]),2),'k')
    % 
    % clf
    % scatter(uvnode(38938,1),uvnode(38938,2),'+r')
    % hold on
    % plot(nodexy(trinodes(38938,[1:end 1]),1),nodexy(trinodes(38938,[1:end 1]),2),'r')
    % scatter(uvnode([38937 39693 38935],1),uvnode([38937 39693 38935],2),'+g')
    % plot(nodexy(trinodes(38937,[1:end 1]),1),nodexy(trinodes(38937,[1:end 1]),2),'g')
    % plot(nodexy(trinodes(39693,[1:end 1]),1),nodexy(trinodes(39693,[1:end 1]),2),'g')
    % plot(nodexy(trinodes(38935,[1:end 1]),1),nodexy(trinodes(38935,[1:end 1]),2),'g')

    % Closed boundary nodes
    %xr=ncread(file,'lon_rho');
    %yr=ncread(file,'lat_rho');
    %tTmp=ncread(file,'temp');
    
    
    %xu=ncread(file,'lon_u');
    %yu=ncread(file,'lat_u');
    %uTmp=ncread(file,'u');
    % 
    % p=pcolor(xr,yr,tTmp(:,:,1,1));
    % p.EdgeColor = 'none';
    %p=pcolor(xr,yr,mask_rho);
    %p.EdgeColor = 'none';
    
%     p=pcolor(xu,yu,uTmp(:,:,1,1));
%     p.EdgeColor = 'none';
    % tTmp(1,1,1,1)
    % tTmp(300,1,1,1)
    % tTmp(379,260,1,1)

    tic
    disp('Identifying nodes that are NOT in mesh, and identifying boundary nodes');
    % Identify whether points are on the closed boundary (some non-NaN
    % neighbours) or not (not in the hydrodynamic model domain)
    %tT=tTmp(:,:,1,1);
    %nanIds=find(isnan(tT));
    maskIds=find(mask_rho==0);
    
    %[rowNaN,colNaN]=find(isnan(tT));
    closedBnodes=[];
    notInMeshNodes=[];

    for i=1:length(maskIds)
        if mod(i,10000)==0
            disp([num2str(i) ' ']);
        end
        % Use convolution to find the neighbours of each NaN element; are they NaN?
        %M=zeros(size(tT,1),size(tT,2));
        M=zeros(size(mask_rho,1),size(mask_rho,2));
        %M(nanIds(i)) = 1;
        M(maskIds(i)) = 1;
        %M(rowNaN(i),colNaN(i)) = 1;
        %neighbourVals=tT(conv2(M,[1,1,1;1,0,1;1,1,1],'same')>0);
        neighbourVals=mask_rho(conv2(M,[1,1,1;1,0,1;1,1,1],'same')>0);
        %if any(~isnan(neighbourVals))
        if any(neighbourVals)
            % closed boundary node case
            %closedBnodes=[closedBnodes nanIds(i)];
            closedBnodes=[closedBnodes maskIds(i)];
        else
            % not in mesh case
            %notInMeshNodes=[notInMeshNodes nanIds(i)];
            notInMeshNodes=[notInMeshNodes maskIds(i)];
        end
    end

    % Open boundary nodes
    obc=[1:n1 (n1*2):n1:(n1*n2) n1*n2-1:-1:(n2-1)*n1+1 (n2-2)*n1+1:-n1:n1+1]; % This is the perimeter
%     hold on
%     scatter(nodexy(obc,1),nodexy(obc,2))
%     scatter(nodexy(closedBnodes,1),nodexy(closedBnodes,2))
    toc

    tic
    disp('Removing components of grid that are not in mesh');
    % Then remove the NaN nodes (from everything)
    disp('- thinning corner nodes and making relabelling table');
    nodexyTmp=nodexy; % Remove the relevant rows
    nodexyTmp(notInMeshNodes,:)=[];
    scatter(nodexyTmp(:,1),nodexyTmp(:,2))
    nodexyList=(1:length(nodexy))';
    nodexyList(notInMeshNodes)=[];
    %nodexyList2=zeros(length(nodexy),1);
    [~,nodexyList2]=ismember(1:length(nodexy),nodexyList); % Make a list of the indices of old nodes in the new reduced list (zero elsewhere)
    nodexyList2=nodexyList2';
    
    disp('- thinning trinodes');
    trinodesTmp=trinodes; % If any of the values is in "notInMesh", remove row. Then reduce values by 1 for each row removed prior to a value in the array
    toRem=zeros(length(trinodes),1);
    for i=1:length(trinodes)
        if any(ismember(trinodes(i,:),notInMeshNodes))
            toRem(i)=1;
        end
    end
    trinodesTmp(toRem==1,:)=[]; % Remove the rows
    disp('- relabelling trinodes');
    for row=1:length(trinodesTmp)
        for col=1:3
            trinodesTmp2(row,col)=nodexyList2(trinodesTmp(row,col));
        end
    end

    disp('- thin neighbouring elements');
    nbeTmp=nbe;
    elementList=(1:length(trinodes))';
    elementList(toRem==1,:)=[];
    nbeTmp(toRem==1,:)=[];
    [~,elementList2]=ismember(1:length(trinodes),elementList); % Make a list of the indices of old elements in the new reduced list (zero elsewhere)
    elementList2=elementList2';
    for row=1:length(nbeTmp)
        for col=1:3
            if nbeTmp(row,col)~=0
                nbeTmp2(row,col)=elementList2(nbeTmp(row,col));
            else
                nbeTmp2(row,col)=0;
            end
        end
    end
    % Elements on the boundary could have a neighbour listed where it
    % should now be zero. Or could they? When putting the values back in,
    % the removed elements should be replaced by zero. Plot to check.
        
    disp('- recreating centroids');
    % e.g. if remove 3,4,5, and a value in trinodes is 6, reduces to 3
    %uvnode % Create again based on the revised nodexy and trinodes
    uvnode=zeros(length(trinodesTmp),2);
    for i=1:length(uvnode)
        uvnode(i,1)=mean(nodexy(trinodesTmp(i,:),1));
        uvnode(i,2)=mean(nodexy(trinodesTmp(i,:),2));
    end
    
    disp('- thinning closed boundary nodes');
    closedBnodesTmp=nodexyList2(closedBnodes)'; % for closed boundary nodes, take on the new node id
    
    disp('- thinning open boundary nodes');
    obcTmp=obc; % Remove values which are NaN
    toRem=zeros(length(obcTmp),1);
    for i=1:length(obcTmp)
        if ismember(obcTmp(i),notInMeshNodes)
            toRem(i)=1;
        end
    end
    obcTmp(toRem==1)=[];
    for i=1:length(obcTmp)
        obcTmp(i)=nodexyList2(obcTmp(i));
    end
    toc
    
    %% Save Mesh
    tic
    disp('Creating Mesh struct');
    Mesh.nodexy=nodexyTmp;
    Mesh.nodexy_deg=nodexyTmp;
    Mesh.uvnode=uvnode;
    Mesh.uvnode_deg=uvnode;
    Mesh.trinodes=trinodesTmp2;
    Mesh.obc=obcTmp;
    Mesh.closedBnodes=closedBnodesTmp;
    Mesh.nbe=nbeTmp2;
    
%     file='D:\hydroOut\NEA_ROMS\3hr\2019\NEATL_SAMS_2019043021.nc';
%     lon_rho=ncread(file,'lon_rho');
%     lat_rho=ncread(file,'lat_rho');
%     depth_rho=ncread(file,'h');
    Fu = scatteredInterpolant(reshape(lon_rho,[],1),reshape(lat_rho,[],1),reshape(depth_rho,[],1));
    Mesh.uvdepth = Fu(Mesh.uvnode(:,1),Mesh.uvnode(:,2));
    Mesh.depth = Fu(Mesh.nodexy(:,1),Mesh.nodexy(:,2));

    %file='D:\hydroOut\NEA_ROMS\1hr\2019\NEATL_2019102300.nc'
    Cs_w=ncread(file,'Cs_w');
    %s_rho=ncread(file,'s_rho');
    % Find boundaries for bottom 11, middle 20, top 10
    Cs_w=Cs_w([1,11,31,37,41]);
    % Get the new layer mid depths
    Mesh.sigvec=Cs_w(1:4)+diff(Cs_w)/2;
    
    toc
    

end

