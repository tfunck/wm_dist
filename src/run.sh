rm -r build/ __pycache__/  /home/t/.local/lib/python3.6/site-packages ; 
python3 setup.py install --user --force
#gdb  -ex=run --args  
python3 BrainDistRun.py --classified "CTRL_C01_pve_classify.mnc"  --surface "CTRL_C01_white_surface_left_81920.obj" --label 3  --matrix  "matrix.csv" --max-threads 4 --subsample "5,4" 
