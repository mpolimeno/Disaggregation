#! /bin/bash

nohup /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nodesktop -nosplash -r "CollectStatsForRotationalFlowIAA; exit" > logs.out &
