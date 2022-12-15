import os, sys, shutil
from datetime import datetime, timedelta

# Read the first command-line argument as the source folder
source_dir = sys.argv[1]

# Read the second command-line argument as the name of the dest folder
dest_dir = sys.argv[2]

# Set start date YYYY/MM/DD
start_date = sys.argv[3]

# Set interval in days
days_interval = int(sys.argv[4])

files = set()
for r,s,f in os.walk(source_dir):
    for file in f:
        files.add(os.path.join(r,file))

print(f"{len(files)} files found in {source_dir}")

range_start = datetime.strptime(start_date, "%Y/%m/%d")
range_end = range_start + timedelta(days=days_interval)

while True:
    tranche_dir = "tranche_{}_{}".format(range_start.strftime("%Y%m%d"), range_end.strftime("%Y%m%d"))
    print(f"tranche dir: {tranche_dir}")
    
    files_in_interval = {f for f in files if os.path.getctime(f) > datetime.timestamp(range_start) and os.path.getctime(f) < datetime.timestamp(range_end)}
    print(f"{len(files_in_interval)} files in this interval")

    for x in files_in_interval:
        subdir = os.path.dirname(x).replace(f"{source_dir}/", '')
        target_dir = os.path.join(dest_dir,tranche_dir,subdir)
        #print(f"file: {x}, subdir: {subdir}, target_dir: {target_dir}")
        os.makedirs(target_dir,exist_ok=True)
        shutil.move(x, target_dir)
    files = files - files_in_interval

    if range_end > datetime.now(): break

    range_start = range_end
    range_end = range_start + timedelta(days=days_interval)