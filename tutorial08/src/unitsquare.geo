// mesh width associated with points
h_corner = 0.0025;
h_top = 0.005;
h_closecorner = 0.008;
h_rest = 0.009;

Point(1) = {0, 0, 0, h_rest};
Point(2) = {1, 0, 0, h_rest};
Point(3) = {1, 0.8, 0, h_closecorner};
Point(4) = {1, 1, 0, h_corner};
Point(5) = {0.8, 1, 0, h_top};
Point(6) = {0.2, 1, 0, h_top};
Point(7) = {0, 1, 0, h_corner};
Point(8) = {0, 0.8, 0, h_closecorner};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,1};

Line Loop(100) = {1,2,3,4,5,6,7,8};  
Plane Surface(200) = {100};  
