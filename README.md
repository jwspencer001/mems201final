# mems201final
To use the code, input your parameters into the variables of the struct defined in “define_truss_problem”, and reference the name of your struct in the truss solver using the line “TR = define_truss_problem”. If you change the name of your struct, be sure to change the name of the reference in the truss solver as well.

The parameters you will need before using this program are:<br/>
-The truss geometry (node positions, members, etc.)<br/>
-Which nodes are fixed in x and or y<br/>
-The forces on each node in x and y<br/>
-The cross sectional area of the members<br/>
-The elastic modulus of the material used for the members

To input the initial positions, put the x position of that node in TR.x0 and the y position in TR.y0. The index of each node should match up in these two vectors. To create the links between nodes (members), put the number of the two nodes that are connected on a new row. The number of rows in this matrix should match the number of members and the number of columns should always be 2. For TR.EA, input the elastic modulus of the material the truss is constructed of multiplied by the cross sectional area of the truss member. For the fixed nodes, using the same format as used to input them initially, enter either 1 or 0, where 1 is a fixed node and 0 is an unfixed node. To apply external forces at the nodes, using the same format as TR.x0/TR.y0, input the x and y force magnitude into TR.loadVectorX and TR.loadVectorY respectively.

How to read the console output: The displacement and nodeForces outputs are in x;y;x;y format (ie. the 2n-1 row is the x value for node n and the 2n row is the y value for node n). The axial forces are in the order that they were entered initially.

Assumptions:<br/>
-Each member of the truss has the same cross sectional area<br/>
-Each member is made of the same material<br/>
-Each joint is essentially treated as a hinge; no fixed connections at the ends of the beams are assumed

Note:<br/>
-Buckling failure is not accounted for<br/>
-Shear failure at the pins is not accounted for<br/>
-Units must be consistent<br/>
-No error checking for the positions and connections is implemented, make sure your inputs are physically possible using the displacement graph.
