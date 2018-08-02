#rm -r build/ __pycache__/  /home/t/.local/lib/python3.6/site-packages ; 
#python3 setup.py install --user --force
#gdb  -ex=run --args  
python3 wm_dist.py --classified "cls.mnc"  --surface "mid.obj" --label 3  --matrix  "matrix.csv" --max-threads 4 --subsample "5,4" 
