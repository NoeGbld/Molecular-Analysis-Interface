import matplotlib
matplotlib.use('QtAgg')

import matplotlib.pyplot as plt
import numpy as np

def plot_absorption_spectrum(mol):
    # 1. Validate incoming data attributes
    if (not hasattr(mol, 'spectre') or mol.spectre is None or
        not hasattr(mol.spectre, 'nu_GHz') or mol.spectre.nu_GHz is None or len(mol.spectre.nu_GHz) == 0 or
        not hasattr(mol.spectre, 'coef') or mol.spectre.coef is None or len(mol.spectre.coef) == 0):
        print("[Plot] No spectral data available")
        return

    freqs = np.array(mol.spectre.nu_GHz)
    coefs = np.array(mol.spectre.coef)

    # 2. Filter data near peak threshold
    if len(coefs) > 0:
        max_coef = np.max(coefs)
        mask = coefs >= (max_coef - 5)
        freqs, coefs = freqs[mask], coefs[mask]

    if len(coefs) == 0:
        print("[Plot] No points to plot after filtering.")
        return

    # Identification de la base de données source (CDMS, JPL ou HITRAN)
    db_name = getattr(mol.spectre, 'databank', 'Unknown')
    db_lower = db_name.lower()
    
    if 'cdms' in db_lower:
        db_label = 'CDMS'
    elif 'jpl' in db_lower:
        db_label = 'JPL'
    elif 'hitran' in db_lower:
        db_label = 'HITRAN'
    else:
        db_label = db_name

    # 3. Handle interactive Spyder GUI canvas generation
    plt.figure(figsize=(10, 4))
    mol_name = getattr(mol, 'name', None) or getattr(mol, 'formula', 'Unknown')
    
    # Remplacement du scatter par un tracé de lignes verticales (plus adapté pour du discret)
    plt.scatter(freqs, coefs, color='red', s=3, marker='o', alpha=0.7, label='Points')
    
    plt.xlabel('Frequency (GHz)')
    
    # Ajustement dynamique de l'unité de l'axe Y selon la base de données
    if 'hitran' in db_lower:
        plt.ylabel(r"Intensity ($\log(\mathrm{cm}^{-1}/(\mathrm{molecule}\cdot\mathrm{cm}^{-2}))$)")
    else:
        plt.ylabel(r"Intensity ($\log_{10}(\mathrm{nm}^2\cdot\mathrm{MHz})$)")
        
    plt.title(f'Discrete Resonance Points of {mol.formula} ({db_label} - Top 5 Decades)')
    plt.grid(True, linestyle='--', alpha=0.5)
    
    # Ajustement des marges pour éviter que les lignes ne touchent les bords
    plt.ylim(np.min(coefs) - 0.5, np.max(coefs) + 0.5)
    plt.tight_layout()
    
    # Render interactive figure modal window
    plt.show()