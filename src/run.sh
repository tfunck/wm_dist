rm -r build/ __pycache__/  /home/t/.local/lib/python3.6/site-packages/BrainDist.cpython-36m-x86_64-linux-gnu.so  ; 
python3 setup.py install --user --force
#gdb  -ex=run --args  
#python3 BrainDistRun.py --classified "CTRL_C01_pve_classify.mnc"  --surface "CTRL_C01_white_surface_left_81920.obj" --label 3  --matrix  "matrix.csv" --max-threads 4 --subsample "5,4" 

gdb  -ex=run --args  python3 surf_dist.py --surface  "mid.obj" --mask "lbl.txt"  --output "surf_dist_test.csv"
#gdb  -ex=run --args  python3 wm_dist.py --classified "cls.mnc" --surface  "mid.obj" -S "lbl.txt"  --matrix "surf_dist_test.csv"
