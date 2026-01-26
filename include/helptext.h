#pragma once

static const char* HELP_TEXT = R"(Usage: ./LIS [OPTIONS]

Simulation of animal interactions in a pen.

Options:
  -h, -H, --help          Show this help message and exit
  -u, -U, --update        Automatically update from github and rebuild 
  -v, -V, --version       Print program version and exit
  
  -i /PATH/TO/INPUT_FILE.ini    Input parameter file
  -o /PATH/TO/OUTPUT_FOLDER/   Output FOLDER
  --seed INTEGER                     Integer random seed (default: time-based)
  --step INTEGER                Number of simulation steps
  --checkpoint INTEGER          Write output every M steps
  --tunnel TRUE/FALSE   TRUE = create tunnel, FALSE = no tunnels
  --mirror TRUE/FALSE   TRUE = mirror, FALSE = no mirror
   


Example:
  BQBN-Reborn.exe -i parameters.ini -o ./experiment1 -seed 123456 --step 100000 --checkpoint 300 --tunnel TRUE --mirror FALSE

Report bugs to: wangzs1995@163.com
)";