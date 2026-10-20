# -*- coding: utf-8 -*-

import os
import sys
import subprocess

def install_dependencies():
    requirements_path = os.path.join(os.path.dirname(__file__), 'requirements.txt')
    
    if not os.path.exists(requirements_path):
        print("[Error] 'requirements.txt' file not found in this directory.")
        return False

    print("==========================================================")
    print(" FIRST TIME SETUP: Installing Required Libraries          ")
    print(" (Isolated User Space - University Security Friendly)      ")
    print("==========================================================")
    
    try:
        print("[1/2] Upgrading local pip package manager...")
        subprocess.run([
            sys.executable, "-m", "pip", "install", 
            "--upgrade", "pip", 
            "--user",
            "--no-warn-script-location"
        ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        print("[2/2] Installing packages (RADIS, Astroquery, Flask, etc.)...")
        print("Please wait, this might take 1 to 2 minutes depending on your network...")
        
        subprocess.run([
            sys.executable, "-m", "pip", "install",
            "--user", 
            "--ignore-installed", 
            "--no-warn-script-location", 
            "-r", requirements_path
        ], check=True)
        
        print("\n[Success] Installation completed successfully!")
        print("You can now close this window and run 'lancer.py' / 'launcher.py'.")
        print("----------------------------------------------------------")
        return True
        
    except Exception as e:
        print(f"\n[Failure] Installation encountered an error: {e}")
        print("Please check your internet connection or university proxy settings.")
        print("----------------------------------------------------------")
        return False

if __name__ == "__main__":
    install_dependencies()
    input("\nPress Enter to exit...")