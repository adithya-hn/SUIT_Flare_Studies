import subprocess
import os

scripts = [
    
   
    ("/media/adithya/Adi_disk4/SUIT_flare_work/case9_nov13/scripts", "aia_full_lc.py"),
    ("/media/adithya/Adi_disk4/SUIT_flare_work/case10_nov13/scripts", "aia_full_lc.py"),
    
]

for folder, script in scripts:
    print(f"Running {script} in {folder}...")
    subprocess.run(["python3", script], cwd=folder, check=True)
    print(f"Finished {script}\n")
    print('--------------------//////\\\\\\---------------------')

'''
 ("/media/adithya/Adi_disk4/SUIT_flare_work/case8_nov01/scripts", "aia_full_lc.py"),
("/media/adithya/Adi_disk4/SUIT_flare_work/case5_jul10/scripts", "aia_full_lc.py"),
    ("/media/adithya/Adi_disk4/SUIT_flare_work/case6_oct09/scripts", "aia_full_lc.py"),
    ("/media/adithya/Adi_disk4/SUIT_flare_work/case7_nov01/scripts", "aia_full_lc.py"),'''