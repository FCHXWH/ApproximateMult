# OPACT: Optimization of Approximate Compressor Tree for Approximate Multiplier 

## Folders

### Src
All .cpp files of this project;
### Header
All .h files of this project;
### jsoncpp-1.9.4
The library **Jsoncpp** is used in this project;
### Verilog
Contains all verilog files (bit_width/OPACT1(2).v, MultiDesignComparison/OPACTDesign1(2,3...).v) synthesized by OPACT in the experimental section of the paper.

## Prerequisities

### Window or Ubuntu
1. Install the Gurobi 9.5;

2. Download the library **Jsoncpp** (just download the folder **jsoncpp-1.9.4** and set the directory **your_directory/jsoncpp-1.9.4/include** during the compilation of this project);

3. Download EDA tool **yosys** (https://yosyshq.net/yosys/download.html) and put the yosys executable file into the same folder of this project's source files;

4. 5. Install Python HDL library **MyHDL**:

    Download this lib in python
    ```bash
    pip install Myhdl
    ```

5. Install the verilog generator:

    Download the source codes of **ApproxMult_Myhdl**:
    ```bash
    git clone https://github.com/FCHXWH/ApproxMULT_MyHDL.git
    ```

6. Install the simulator (option):
    
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
    


