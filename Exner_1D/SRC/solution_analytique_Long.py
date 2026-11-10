import numpy as np
from scipy.integrate import simps
from scipy.optimize import fsolve

# --- Paramètres du problème ---
a, Q, np_poro, s = 0.001, 10.0, 0.4, 0.0
h0, beta, xc = 6.0, 0.01, 150.0
T = 6000.0  # Modifier ici pour 2000.0

# Seuil de détection du choc (en mètres)
# Si le repli mathématique est plus petit que ça, on reste en solution forte.
overlap_threshold = 0.05

def C_phys(zb):
    return (3 * a * Q**3) / ((1 - np_poro) * (s - zb)**4)

# 1. Génération de la solution "Forte"
x0 = np.linspace(0, 300, 10000)
zb0 = -h0 + 2.0 * np.exp(-beta * (x0 - xc)**2)
x_strong = x0 + C_phys(zb0) * T

# 2. Détection intelligente du choc
# On calcule le repli maximum (valeur la plus négative du saut de x)
diff_x = np.diff(x_strong)
max_overlap = -np.min(diff_x) if np.any(diff_x < 0) else 0.0

if max_overlap < overlap_threshold:
    # --- SOLUTION FORTE ---
    print(f"T={T}s : Repli ({max_overlap:.4f}m) sous le seuil. Solution forte utilisée.")
    x_res = x_strong
    z_res = zb0
else:
    # --- SOLUTION FAIBLE (Choc) ---
    print(f"T={T}s : Choc détecté (Repli={max_overlap:.2f}m). Calcul de la solution faible.")
    aire_ref = simps(zb0 + h0, x0)

    def find_shock(xs):
        # On ne regarde que la face arrière pour définir la forme avant le choc
        mask = (x0 <= xc) & (x_strong <= xs)
        if not np.any(mask): return aire_ref
        x_val, z_val = x_strong[mask], zb0[mask]
        current_aire = simps(z_val + h0, x_val)
        # On ajoute le rectangle entre le dernier point transporté et la verticale du choc
        current_aire += (z_val[-1] + h0) * (xs - x_val[-1])
        return current_aire - aire_ref

    # Recherche de la position du choc
    x_shock = fsolve(find_shock, xc + 10.0)[0]
    
    # Construction finale
    mask_back = (x0 <= xc) & (x_strong <= x_shock)
    x_final = x_strong[mask_back]
    z_final = zb0[mask_back]
    
    x_res = np.concatenate([x_final, [x_shock, x_shock, 300.0]])
    z_res = np.concatenate([z_final, [z_final[-1], -h0, -h0]])

# 3. Exportation
output_path = f"./SORTIE/LONG/ana_faible_{int(T)}.dat"
np.savetxt(output_path, np.column_stack((x_res, z_res)), header="x zb")
print(f"Fichier généré : {output_path}")