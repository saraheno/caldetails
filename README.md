#Single Bar
GEANT4-based simulation of a single bar for DualReadout purpose.

## Installing

Currently this only works on UMD cluster. Might need some modification if using on other clusters. 
```bash
source g4env.sh
cmake -DGeant4_DIR=/cvmfs/geant4.cern.ch/geant4/10.5/x86_64-slc6-gcc63-opt/lib64/GEANT4-10.5.0
cmake --build .
```

## Single event running
For the interactive runing:
Check the geometry through a detector viewing window
```bash
./CEPC_CaloTiming -c template.cfg -u Xm  
```
Run the events defined in run.mac
```bash
./CEPC_CaloTiming -c template.cfg -m run.mac -o test
# -c config file
# -m particle source
# -o output file
```

## Bulk running
use the generate_temp.sh in `./template`
```bash
cd template 
source generate_temp.sh
```




