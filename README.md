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
    
    


