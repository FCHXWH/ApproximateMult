# OPACT: Optimization of Approximate Compressor Tree for Approximate Multiplier 

## Folders

### Src
All .cpp files of this project;
### Header
All .h files of this project;
### jsoncpp-1.9.4
Download jsoncpp library (https://github.com/open-source-parsers/jsoncpp)
### Verilog
Contains all verilog files (bit_width/OPACT1(2).v, MultiDesignComparison/OPACTDesign1(2,3...).v) synthesized by OPACT in the experimental section of the paper.

## Prerequisities

### Window or Ubuntu
1. Install the Gurobi 9.5;

2. Download EDA tool **yosys** (https://yosyshq.net/yosys/download.html) and put the yosys executable file into the same folder of this project's source files;

3. Install Python HDL library **MyHDL**:

    Download this lib in python
    ```bash
    pip install Myhdl
    ```

4. Install the verilog generator:

    Download the source codes of **ApproxMult_Myhdl**:
    ```bash
    git clone https://github.com/FCHXWH/ApproxMULT_MyHDL.git
    ```

5. Install the simulator (option):
    
    Download the source codes of **simulator**
    ```bash
    git clone https://github.com/changmg/simulator.git
    ```
    
    Compile
    ```bash
    cd simulator/
    mkdir build
    cd build
    cmake ..
    make -j16
    cd ..
    ```
    How to use this simulator (simulation.out -h):  
   simulation.out -a exact_multiplier.blif -b approximate_multiplier.blif -f error_report_file_path
    

## Compilation

### Include directory
1. your_directory/jsoncpp-1.9.4/include;
2. your_directory/Gurobi/win64/include.

### Library directory
1. your_directory/Gurobi/win64/lib;
2. your_directory/jsoncpp-1.9.4/jsonlib.

### Linker
1. gurobi_c++mdd2019.lib;
2. gurobi91.lib;
3. jsoncpp.lib.

### Configuration in the code
1. _Opt_Mode_: 0=>AreaOpt; 1=>ErrOpt; 2=>CoOpt;
2. _bounds_: set the constraint of different metrics under different modes (for example, in AreaOpt mode, set the bound of error);
3. _path_: set the directory of the downloaded verilog generator.
