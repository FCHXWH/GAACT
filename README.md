# GAACT: Genetic Algorithm for Approximate Compressor Tree 
## Prerequisities

### Ubuntu
1. Install Berkely ABC;

2. Install the lastest linux Gurobi 9.5;

3. Install the library **Jsoncpp**:

    Download the source codes of **Jsoncpp**
    ```bash
    git clone https://github.com/open-source-parsers/jsoncpp/tree/0.y.z;
    ```
    Compile
    ```bash
    unzip jsoncpp-0.y.z.zip
    cd jsoncpp-0.y.z
    mkdir -p build/debug
    cd build/debug
    cmake -DCMAKE_BUILD_TYPE=debug -DBUILD_STATIC_LIBS=ON -DBUILD_SHARED_LIBS=ON -DARCHIVE_INSTALL_DIR=. -G "Unix Makefiles" ../..
    make
    ```

4. Install the simulator:
    
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
    
5. Install Python HDL library **MyHDL**:

    Download this lib in python
    ```bash
    pip install Myhdl
    ```

6. Install the verilog generator:

    Download the source codes of **ApproxMult_Myhdl**:
    ```bash
    git clone https://github.com/FCHXWH/ApproxMULT_MyHDL.git
    ```
## What we have done
- Encoding a connection order as a gene through **Permutation Encoding**;
- Define the self-designed crossover, mutation operators, which are easy to be parallelized;
- Learn from the implementation of **priority cuts** in **Berkely ABC**:
    1. Implement the incremental sorting function to accelerate the selection operator;
    2. Assign a signature for each gene;
    3. Set up a hash table for each population to prune the repeated genes to optimize the memory management;
- Set up a global hash table to record the times of each gene, which is the key of UCB algorithm.

## TBD
- [x] Further optimize the memory management of genes;
- [x] Combine the EE (Exploration-Exploitation) trade-off with the **Selection** operator in GA;
- [ ] Implement NSGA-II algorithm to co-optimize the delay and error of approximate CT;
- [ ] Further parallelize GA operations: crossover, mutation, etc; 
- [ ] Further improve the accuracy of error estimation.
    


