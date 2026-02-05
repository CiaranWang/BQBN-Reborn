#pragma once

static const char* HELP_TEXT = R"(Usage: BQBN-Reborn.exe [OPTIONS]

Simulation of LBM.

Options:
  -h, -H, --help          Show this help message and exit
  -u, -U, --update        Automatically update BQBN-Reborn.exe from Github 
  -v, -V, --version       Print program version and exit
  
  -i /PATH/TO/INPUT_FILE.ini    Input parameter file
  -o /PATH/TO/OUTPUT_FOLDER/   Output FOLDER
  -n NAME                       experiment name, this will be output file name prefix (e.g. NAME.info)
  --seed INTEGER                Integer random seed (default: time-based)
  --step INTEGER                Number of simulation steps (default: 1000)
  --checkpoint INTEGER          Write output every M steps (default: 200)
  
  --tunnel [VALUE]               Enable or disable tunnel mode (default: false)
  --mirror [VALUE]               Enable or disable mirror mode (default: false)

    VALUE interpretation:
      1, true, TRUE, or no value      → TRUE
      0, false, FALSE                 → FALSE
      (flag not provided)             → FALSE (default)
   


Example:
  BQBN-Reborn.exe -i parameters.ini -o ./exp1 -seed 123 --step 10000 --checkpoint 300 --tunnel --mirror
  BQBN-Reborn.exe -i parameters.ini -o ./exp1 -seed 123 --step 10000 --checkpoint 300 --tunnel 1 --mirror true
  BQBN-Reborn.exe -i parameters.ini -o ./exp1 -seed 123 --step 10000 --checkpoint 300
Report bugs to: wangzs1995@163.com
)";