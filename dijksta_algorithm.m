function [dist,path] = dijksta_algorithm(verteces,edges,start_vertex,destination_vertex)
% 

if (nargin < 3) % SETUP
%     NODES CONSTRUCTING
    num_verteces = 83; 
    coordinates=[1 1;1 2;1 3;1 4;1 5;1 6;1 7;1 8;1 9;1 10;
           2 1;2 7;2 8;2 9;2 10;
           3 1;3 2;3 3;3 4;3 5;3 6;3 7;3 8;3 9;3 10;
           4 1;4 2;4 3;4 4;4 5;4 6;4 7;4 8;4 9;4 10;
           5 1;5 2;5 3;5 4;5 6;5 7;5 8;5 9;5 10;
           6 1;6 2;6 3;6 4;6 6;6 7;6 8;6 9;6 10;
           7 1;7 2;7 3;7 4;7 5;7 6;7 7;7 8;7 9;7 10;
           8 1;8 2;8 3;8 4;8 5;8 6;8 7;8 8;8 9;8 10;
           9 1;9 2;9 3;9 4;9 5;9 6;9 7;9 8;9 9;9 10];
          
    max_edge_length = 1.4;
    ids = (1:num_verteces)';
    verteces = [ids coordinates]; % create random verteces
    h = figure; plot(verteces(:,2),verteces(:,3),'k.') % plot the verteces
    text(verteces(num_verteces,2),verteces(num_verteces,3),...
        [' ' num2str(ids(num_verteces))],'Color','b','FontWeight','b')
    hold on
%     EDGES CONSTRUCTING
    num_edges = 0; edges = zeros(num_verteces*(num_verteces-1)/2,3);
    for i = 1:num_verteces-1 % create edges between some of the verteces
        text(verteces(i,2),verteces(i,3),[' ' num2str(ids(i))],'Color','b','FontWeight','b')
        for j = i+1:num_verteces
            d = sqrt(sum((verteces(i,2:3) - verteces(j,2:3)).^2));
                if d < max_edge_length
                plot([verteces(i,2) verteces(j,2)],[verteces(i,3) verteces(j,3)],'k.-')
                % add this link to the edges list
                num_edges = num_edges + 1;
                edges(num_edges,:) = [num_edges verteces(i,1) verteces(j,1)];
            end
        end
    end
    edges(num_edges+1:num_verteces*(num_verteces-1)/2,:) = [];
    axis([0 11 0 11])
    
    % Calculate Shortest Path Using Dijkstra's Algorithm
    start_vertex = input('Enter start vertex:\n'); 
    disp(['start vertex = ' num2str(start_vertex)]);
    destination_vertex = input('Enter end vertex:\n'); 
    disp(['destination_vertex = ' num2str(destination_vertex)]);
    [distance,path] = dijksta_algorithm(verteces,edges,start_vertex,destination_vertex);
    disp(['distance = ' num2str(distance)]); 
    disp(['path = [' num2str(path) ']']);
    % If a Shortest Path exists,Plot it on the Map.
    figure(h)
    for k = 2:length(path)
        m = find(verteces(:,1) == path(k-1));
        n = find(verteces(:,1) == path(k));
        plot([verteces(m,2) verteces(n,2)],[verteces(m,3) verteces(n,3)],'ro-','LineWidth',2);
    end
    title(['Shortest Distance from ' num2str(start_vertex) ' to ' ...
        num2str(destination_vertex) ' = ' num2str(distance)])
    hold off
    
else %--------------------------------------------------------------------------
    % MAIN FUNCTION - DIJKSTRA'S ALGORITHM
    
    % initializations
    vertex_ids = verteces(:,1);
    [num_map_pts,cols] = size(verteces);
    table = sparse(num_map_pts,2);%sparse creates adjecency matrix
    shortest_distance = Inf(num_map_pts,1);%initializing 
    settled = zeros(num_map_pts,1);%initializing 
    path = num2cell(NaN(num_map_pts,1));
    col = 2;
    pidx = find(start_vertex == vertex_ids);
    shortest_distance(pidx) = 0;
    table(pidx,col) = 0;%creating table with dimensions pidx*columns(=2)
    settled(pidx) = 1;
%     we checked start point
    path(pidx) = {start_vertex};
    if (nargin < 4) % compute shortest path for all verteces
        queing = 'sum(~settled) > 0';
    else % terminate algorithm early
        queing = 'settled(zz) == 0';
        zz = find(destination_vertex == vertex_ids);
    end
    while eval(queing)
        % update the table
        table(:,col-1) = table(:,col);
        table(pidx,col) = 0;
        % find neighboring verteces in the edges list
        neighbour_ids = [edges(vertex_ids(pidx) == edges(:,2),3);
            edges(vertex_ids(pidx) == edges(:,3),2)];
        % calculate the distances to the neighboring verteces and keep track of the paths
        for k = 1:length(neighbour_ids)
            cidx = find(neighbour_ids(k) == vertex_ids);
            if not(settled(cidx))%for zero elements in settled
                d = sqrt(sum((verteces(pidx,2:cols) - verteces(cidx,2:cols)).^2));
                if (table(cidx,col-1) == 0) || ...
                        (table(cidx,col-1) > (table(pidx,col-1) + d))
                    table(cidx,col) = table(pidx,col-1) + d;
                    tmp_path = path(pidx);
                    path(cidx) = {[tmp_path{1} neighbour_ids(k)]};
                else
                    table(cidx,col) = table(cidx,col-1);
                end
            end
        end
        % find the minimum non-zero value in the table and save it
        nidx = find(table(:,col));
        ndx = find(table(nidx,col) == min(table(nidx,col)));
        if isempty(ndx)
            break
        else
            pidx = nidx(ndx(1));
            shortest_distance(pidx) = table(pidx,col);
            settled(pidx) = 1;
        end
    end
    if (nargin < 4) % return the distance and path arrays for all of the verteces
        dist = shortest_distance';
        path = path';
    else % return the distance and path for the ending node
        dist = shortest_distance(zz);
        path = path(zz);
        path = path{1};
    end
end
