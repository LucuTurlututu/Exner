// On définit l'algorithme de maillage 2D sur 5 (Delaunay)
// Les valeurs possibles sont : 1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 8=DelQuad
Mesh.Algorithm = 5; 

// On définit la taille caractéristique des mailles (lc)
// Pour avoir environ la même résolution que ton "N = 160" sur un domaine de 7x7 :
// lc = Longueur / N = 7 / 160 ≈ 0.04375
lc = 1/160;

// Points avec la taille de maille cible lc
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};

Physical Surface("Domaine", 1) = {1};
