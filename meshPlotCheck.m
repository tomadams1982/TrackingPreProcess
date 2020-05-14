function meshPlotCheck(Mesh)
%% Plot things up to check
    disp('Plotting');
    clf
    plotMeshPDens(Mesh,'os',0,'plotEdges',1)
    hold on

    try
        scatter(Mesh.nodexy(Mesh.closedBnodes,1),Mesh.nodexy(Mesh.closedBnodes,2),'b')
    catch
        try
            scatter(Mesh.nodexy(Mesh.boundaryNodesAll,1),Mesh.nodexy(Mesh.boundaryNodesAll,2),'b')
        catch
            warning('Fields closedBnodes or boundaryNodesAll not found');
        end
    end
    try
        scatter(Mesh.nodexy(Mesh.obc,1),Mesh.nodexy(Mesh.obc,2),'r')
    catch
        warning('Field obc was not present');
    end

    box on
    %print('-painters','-dpng','-r600','NEATL_SAMS_domainCloseUp2.png')
    % plot bElems
    bElem=zeros(length(Mesh.uvnode),1);
    for row=1:length(Mesh.uvnode)
        if any(Mesh.nbe(row,:)==0)
            bElem(row)=1;
        end
    end
    scatter(Mesh.uvnode(bElem==1,1),Mesh.uvnode(bElem==1,2),'+g');
    
%     scatter(lon_rho(200,200),lat_rho(200,200),'r')
%     scatter(lon_rho(800,200),lat_rho(800,200),'r')
%     scatter(lon_rho(200,750),lat_rho(200,750),'r')
%     scatter(lon_rho(800,750),lat_rho(800,750),'r')
end