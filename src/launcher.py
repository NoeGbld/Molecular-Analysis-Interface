# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import webbrowser
import threading
import time
import socket

def open_browser():
    # Wait until Flask answers (first launch downloads Chromium, this can take minutes)
    for _ in range(1500):  # ~5 min max
        try:
            with socket.create_connection(("127.0.0.1", 5000), timeout=0.2):
                break
        except OSError:
            time.sleep(0.2)
    print("[System] Automatically opening the web browser...")
    webbrowser.open_new("http://127.0.0.1:5000/")

if __name__ == "__main__":
    print("[System] Starting the application...")
    
    # 1. Trigger the browser launch routine in a separate background thread
    threading.Thread(target=open_browser, daemon=True).start()
    
    # 2. Path to the Flask server
    app_path = os.path.join(os.path.dirname(__file__), 'app.py')
    
    if not os.path.exists(app_path):
        print(f"[Error] Target file '{app_path}' is missing.")
        input("Press Enter to exit...")
        sys.exit(1)
        
    # 3. Use Popen instead of run to monitor the background process
    try:
        process = subprocess.Popen([sys.executable, app_path],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           text=True, bufsize=1)
        threading.Thread(target=lambda: [print(l, end="") for l in process.stdout],
                 daemon=True).start()
        
        # Active monitoring loop
        while True:
            # poll() returns None if the process is still running
            if process.poll() is not None:
                print("\n[System] Shutdown detected from web interface.")
                break
            time.sleep(0.5)  # Check every 500ms to keep CPU usage at 0%
            
    except KeyboardInterrupt:
        # Fallback if someone uses Ctrl+C anyway
        print("\n[System] Ctrl+C detected. Forcing shutdown...")
        if process and process.poll() is None:
            process.terminate()
            
    print("[System] Application closed successfully. Safe to close this window.")
    time.sleep(1.0)
    sys.exit(0)