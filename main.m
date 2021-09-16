%% Section 1  - THE Stiffness Matrix

% load file with parameters for the node position, fixed nodes, connections, forces:
    TR = define_truss_problem_final1;
    
% create initial stiffness matrix of proper size filled with zeros
    stiffness = ...
        zeros(length(TR.x0)+length(TR.y0),length(TR.x0)+length(TR.y0));

% Begins with the creation of the local stiffness matrix and 
% transformation of the local stiffness matrix to global coordinates. 
% Proceeds to insert the four portions of the local stiffness matrix in
% global coordinates to the global stiffness matrix.
for i = 1:TR.numEL
    % reference the nodes connected to the selected member.
    node1=TR.links(i,1);
    node2=TR.links(i,2);
    % x and y lengths of the member/bar
    delta_x=TR.x0(node2)-TR.x0(node1);
    delta_y=TR.y0(node2)-TR.y0(node1);
    L = sqrt( (delta_x).^2 + (delta_y).^2 ); % length of the bar 
    c = (delta_x) / L; % cosine of bar angle 
    s = (delta_y) / L; % sine of bar angle 
    % Local stiffness matrix in global coordinates. c and s correspond to
    % the x-proportion and y-proportion of the bar (which are the cosine
    % and sine of the bar angle). The various combinations of c and s in
    % the matrix correspond to their transformations to global coordinates
    % from the local stiffness matrix in local coordinates.
    K = (TR.EA(1)/L) * [ c.^2, c*s, -c.^2, -c*s ;...
                         c*s, s.^2, -c*s, -s.^2 ;...
                        -c.^2, -c*s, c.^2,  c*s ;...
                        -c*s, -s.^2, c*s,  s.^2 ];
    % The following details the superposition of the four quadrants of the
    % local stiffness matrix in global coordinates into the global
    % stiffness matrix "stiffness". The column and row of the global
    % stiffness matrix correspond to the "X" row of the node the local
    % matrix's quadrant must insert into.
        % K top left
    stiffness(2*node1-1:2*node1,2*node1-1:2*node1) = ...
        stiffness(2*node1-1:2*node1,2*node1-1:2*node1) + K(1:2,1:2);
        % K top right
    stiffness(2*node1-1:2*node1,2*node2-1:2*node2) = ...
        stiffness(2*node1-1:2*node1,2*node2-1:2*node2) + K(1:2,3:4);
        % K bottom left
    stiffness(2*node2-1:2*node2,2*node1-1:2*node1) = ...
        stiffness(2*node2-1:2*node2,2*node1-1:2*node1) + K(3:4,1:2);
        % K bottom right
    stiffness(2*node2-1:2*node2,2*node2-1:2*node2) = ...
        stiffness(2*node2-1:2*node2,2*node2-1:2*node2) + K(3:4,3:4);
end

%% Section 2 - Solving

% Copy the global stiffness matrix
kGR = stiffness;

% Create a load vector that is the length of the given X and Y load vectors
%concatenated and then fill it with alternating X and Y values (i.e. x1,y1,x2,y2,etc.)
loadVector = zeros(length(TR.loadVectorX) + length(TR.loadVectorY),1);
loadVector(1:2:end-1) = TR.loadVectorX;
loadVector(2:2:end) = TR.loadVectorY;

% Create a vector of the fixed nodes in the same format as the load vector and the global
%stiffness matrix
fixedNodes = zeros(length(TR.fixedNodesX) + length(TR.fixedNodesY),1);
fixedNodes(1:2:end-1) = TR.fixedNodesX;
fixedNodes(2:2:end) = TR.fixedNodesY;

% Use matlab find in order to find the indices of the fixed nodes
fixedIdx = find(fixedNodes);


% Now eliminate the rows and columns associated with the fixed nodes
kGR(fixedIdx,:) = [];
kGR(:,fixedIdx) = [];

% Copy the load vector
copyLoad = loadVector;

% Eliminate the rows associated with the fixed nodes in the coppied load
%vector similar to the KGR matrix
copyLoad(fixedIdx) = [];

% Solve for the displacements of the unfixed nodes
u_notFixed = kGR\copyLoad;

% This whole section makes a u vector with the displaced and non-displaced points
% Create vector utotal with length that is the number of the nodes times 2 (same as loadVector)
u_total = zeros(length(loadVector),1);
% Add in all of the nodes with a non-zero displacement using a logical array
u_total(find(~fixedNodes)) = u_notFixed;
displacements=u_total

% Now the node reaction force vector for the entire matrix can be found
nodeForces = stiffness * u_total;

% Create a vector of the initial positions in xyxy format
initialPositions = zeros(length(TR.x0) + length(TR.y0),1);
initialPositions(1:2:end-1) = TR.x0;
initialPositions(2:2:end) = TR.y0;

% Now just add the two vectors together to get the final positions
finalPositions = initialPositions + u_total;

%% Section 3 - Ploting

% Plot the initial node positions using the given vectors
plot(TR.x0, TR.y0,'ko','MarkerSize',6);

% Plot the truss elements (the tricky part)
% Plotting of the un-deformed truss members
hold on
for i = 1:TR.numEL
    node1=TR.links(i,1);
node2=TR.links(i,2);
xPos1 = TR.x0(node1);
yPos1 = TR.y0(node1);
xPos2 = TR.x0(node2);
yPos2 = TR.y0(node2);
plot([xPos1 xPos2],[yPos1 yPos2],'k-');
end

%Creating the displaced nodes with an exageration factor
%This changes how much the deformation is multiplied by
exageration_factor = 100;
%This is a copy of the code in the solver section that creates the final
%vector of displaced nodes, but this time multiplied by the exageration
%factor
u_total_exagerated = zeros(length(loadVector),1);
u_total_exagerated(find(~fixedNodes)) = u_notFixed*exageration_factor;

finalPositions_exagerated = initialPositions + u_total_exagerated;

%Unpack the final positions to get back into separate x and y vectors
finalPosX = finalPositions_exagerated(1:2:end-1);
finalPosY = finalPositions_exagerated(2:2:end);

%Plot the final node positions using the found vectors
plot(finalPosX, finalPosY,'ro','MarkerSize',6);

%Plotting of the deformed truss members
for i = 1:TR.numEL
    node1=TR.links(i,1);
node2=TR.links(i,2);
xPos1 = finalPosX(node1);
yPos1 = finalPosY(node1);
xPos2 = finalPosX(node2);
yPos2 = finalPosY(node2);
plot([xPos1 xPos2],[yPos1 yPos2],'r--');
end

%This section creates arrows for the input forces using the matlab quiver
%comand
%setting arrow size through a:
a = 2;

%This loop makes a quiver with a normalized length (controlled by a above) at the position of the
%input forces
for i=1:length(TR.x0)
    quiver(TR.x0(i),TR.y0(i),TR.loadVectorX(i)/abs(TR.loadVectorX(i))*a,0,'k','LineWidth',1.5)
end
for i=1:length(TR.x0)
    quiver(TR.x0(i),TR.y0(i),0,TR.loadVectorY(i)/abs(TR.loadVectorY(i))*a,'k','LineWidth',1.5)
end

%plotting fixed points
for i=1:length(TR.x0)
   if TR.fixedNodesX(i) == 1
       plot(TR.x0(i),TR.y0(i),'.k','MarkerSize',20)
   end
end
hold off

%% Section 4 - Command Window Output

% axialForce determines the axial force in each member by using the change 
% in length, elastic modulus, area, and inital length. The inital lengths 
% and final lengths are taken as the 2-norm of the x and y distances 
% between the two nodes, given by the finalPositions, TR.x0, and TR.y0
% vectors.
axialForce = zeros(TR.numEL,1);
% iterate through each of the members
for i = 1:TR.numEL
    % identify the nodes in contact with the member
	node1 = TR.links(i,1);
	node2 = TR.links(i,2);
    % x and y distances for inital length
    delta_x0 = TR.x0(node2)-TR.x0(node1);
    delta_y0 = TR.y0(node2)-TR.y0(node1);
    % initalLength is the 2-norm of the x and y distances
    initalLength = sqrt( (delta_x0).^2 + (delta_y0).^2 ); 
    % repeat process for the final distances, using the concatenated
    % finalPositions vector
    finalX1 = finalPositions(node1.*2-1);
    finalY1 = finalPositions(node1.*2);
    finalX2 = finalPositions(node2.*2-1);
    finalY2 = finalPositions(node2.*2);
    delta_x = finalX2-finalX1;
    delta_y = finalY2-finalY1;
    finalLength = sqrt(delta_x.^2 + delta_y.^2);
    % calculate the change in length
    delta = finalLength - initalLength;
    % find axialForce of the member using delta
    axialForce(i) = TR.EA(1).*delta./initalLength;
end

%Display the axial forces and the forces at the nodes
axialForce
nodeForces
