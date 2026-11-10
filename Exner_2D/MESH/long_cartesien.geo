//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {400, 0, 0, 1.0};
//+
Point(3) = {400, 1, 0, 1.0};  // Modification de la hauteur Y de 400 à 10
//+
Point(4) = {0, 1, 0, 1.0};   // Modification de la hauteur Y de 400 à 10
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};

// Discrétisation indépendante X et Y
Nx = 300;
Ny = 1;

Transfinite Curve {1, 3} = Nx + 1; // 300 mailles le long de X (lignes horizontales)
Transfinite Curve {2, 4} = Ny + 1; // 10 mailles le long de Y (lignes verticales)

Transfinite Surface {1} = {1, 2, 3, 4};

Recombine Surface {1};

Physical Surface("Domaine", 1) = {1};