 define_truss_problem_final2
% 
% Define the parameters of the truss problem to solve:
% TR = define_truss_problem_final2
% set parameters within this function
% 
% Input:
%     set input parameters within the file
% Output:
%     TR - structure with parameters that define the truss problem
%

% INSTUCTIONS
% To use the code, simply write "TR = define_truss_problem_final2"
% at the beginning of MEMS_201_Final_Project.

% Note: units are not defined below. You will need to keep your force,
% distance, and stiffness units consistent.

% initial position of x,y points
TR.x0 = [0, 3, 6, 0, 3]';
TR.y0 = [0, 0, 0, 3, 3]';

% connections between nodes. 
%    - each row is a truss element
%    - each entry in the row is a node number (the truss element ends)
TR.links = [1 2; 2 3; 4 5; 2 4; 2 5; 5 3]; 

% per element axial rigidity, should be column vector same length as links.
% This is the Young's Modulus E times area A, units are Newtons or lbf.
% Axial stiffness of truss element is then EA/L
TR.EA = 10*10^6*0.1*[1, 1]';

% Fix nodes here
% enter 1 if DOF is fixed, 0 if it is not
TR.fixedNodesX = [1, 0, 0, 1, 0]'; 
TR.fixedNodesY = [1, 0, 0, 1, 0]'; 

% Apply loads here (Newtons)
% enter 0 for no load, or value for load
TR.loadVectorX = [0, 0, 0, 0, 0]'; 
TR.loadVectorY = [0, -200, -200, 0, 0]';

% Apply initial displacements here (mm)
% enter 0 for no prescribed displacement, or value for displacement
% YOU MAY NOT USE THIS, CAN IGNORE
TR.u0VectorX = [0,0,0,0,0]'; 
TR.u0VectorY = [0,0,0,0,0]';

% get the number of nodes and number of truss elements
% NO USER INPUT BELOW
TR.numNodes = length(TR.x0);
TR.numEL = size(TR.links,1);
