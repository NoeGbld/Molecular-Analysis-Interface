import io
import base64
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def plot_absorption_spectrum(mol):
    # Safe check for spectral data availability
    if (not hasattr(mol, 'spectre') or mol.spectre is None or
        not hasattr(mol.spectre, 'nu_GHz') or mol.spectre.nu_GHz is None or len(mol.spectre.nu_GHz) == 0 or
        not hasattr(mol.spectre, 'coef') or mol.spectre.coef is None or len(mol.spectre.coef) == 0):
        print("[Plot] No spectral data available")
        return None
    
    spec_type = getattr(mol.spectre, 'type', 'discrete')
    freqs = np.array(mol.spectre.nu_GHz)
    coefs = np.array(mol.spectre.coef)
    
    # Database identification
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

    plt.figure(figsize=(10, 4))
    mol_name = getattr(mol, 'name', None) or getattr(mol, 'formula', 'Unknown')
    
    # Case 1: Continuous spectrum (RADIS)
    if spec_type == 'continuous':
        plt.plot(freqs, coefs, color='darkblue', linewidth=1.2, label=f'RADIS Continuum ({db_label})')
        plt.ylabel(r"Cross-section ($\mathrm{cm}^2/\mathrm{molec}/\mathrm{cm}^{-1}$)")
            
        plt.title(f'Continuous Absorption Spectrum of {mol.formula} (RADIS - {db_label})')

    # Case 2: Discrete points (Astroquery)
    else:
        # Keep data within the top 5 orders of magnitude
        if len(coefs) > 0:
            max_coef = np.max(coefs)
            mask = coefs >= (max_coef - 5)
            freqs = freqs[mask]
            coefs = coefs[mask]

        # Prevent crash if no points remain after filtering
        if len(coefs) == 0:
            print("[Plot] No discrete points left after filtering top 5 orders of magnitude")
            plt.close()
            return None

        plt.scatter(freqs, coefs, color='red', s=3, marker='o', alpha=0.7, edgecolors='none', label='Resonance Points')
        
        # Ajustement de l'unité de l'axe Y si c'est HITRAN
        if 'hitran' in db_lower:
            plt.ylabel(r"Intensity ($\log(\mathrm{cm}^{-1}/(\mathrm{molecule}\cdot\mathrm{cm}^{-2}))$)")
        else:
            plt.ylabel(r"Intensity ($\log_{10}(\mathrm{nm}^2\cdot\mathrm{MHz})$)")
        
        title_str = f'Discrete Resonance Points of {mol.formula} ({db_label} - Top 5 Decades)'
        plt.title(title_str)
        
        # Dynamic margins based on filtered points
        plt.ylim(np.min(coefs) - 0.5, np.max(coefs) + 0.5)
    
    # Common plot settings
    plt.xlabel('Frequency (GHz)')
    plt.grid(True, linestyle='--', alpha=0.5, which='both' if spec_type == 'continuous' else 'major')
    plt.legend(loc='upper right')
    plt.tight_layout()
    
    # Save to memory buffer
    img = io.BytesIO()
    plt.savefig(img, format='png', dpi=150)
    img.seek(0)
    plt.close()
    
    # Base64 encoding for HTML delivery
    plot_url = base64.b64encode(img.getvalue()).decode('utf-8')
    return f"data:image/png;base64,{plot_url}"