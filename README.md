# GAACT: Genetic Algorithm for Approximate Compressor Tree 
## Prerequisities
### Ubuntu
1. Install the lastest linux Gurobi 9.5:
2. Install the library **Jsoncpp**:
(i). 
```bash
git clone https://github.com/open-source-parsers/jsoncpp/tree/0.y.z;
```
(ii). 
```bash
unzip jsoncpp-0.y.z.zip
cd jsoncpp-0.y.z
```
(iii).
```bash
mkdir -p build/debug
cd build/debug
cmake -DCMAKE_BUILD_TYPE=debug -DBUILD_STATIC_LIBS=ON -DBUILD_SHARED_LIBS=ON -DARCHIVE_INSTALL_DIR=. -G "Unix Makefiles" ../..
make
```



