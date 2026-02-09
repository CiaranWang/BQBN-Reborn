#pragma once

static const char* HELP_TEXT = R"(Usage: BQBN-Reborn.exe [OPTIONS]

Simulation of LBM.

Options:
  -h, -H, --help          Show this help message and exit
  -u, -U, --update        Automatically update BQBN-Reborn.exe from GitHub 
  -v, -V, --version       Print program version and exit
  
  -i /PATH/TO/INPUT_FILE.ini    Input parameter file (REQUIRED)
  -o /PATH/TO/OUTPUT_FOLDER/   Output FOLDER (default: input file folder)
  -n NAME                       experiment name (REQUIRED), this will be output file name prefix. e.g. NAME.info
  --seed INTEGER                Integer random seed (default: time-based)
  --step INTEGER                Number of simulation steps (default: 1000)
  --checkpoint INTEGER          Write output every M steps (default: 200)
  
  --tunnel N               radius of tunnel, tunnel width = 2N+1 (default: no tunnel)
  --mirror [VALUE] Enable mirror mode (VALUE: 1/true/TRUE = on, 0/false/FALSE = off; if omitted => on)   


Example:
  BQBN-Reborn.exe -i parameters.ini -o ./results -n exp1  --seed 123 --step 10000 --checkpoint 300 --tunnel 2 --mirror
  BQBN-Reborn.exe -i parameters.ini -o ./results -n exp1  --seed 123 --step 10000 --checkpoint 300 --mirror true
  BQBN-Reborn.exe -i parameters.ini -n exp1

Report bugs to: wangzs1995@163.com
)";