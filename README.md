# Rainbowfish Color Reuse

This is a project for CIS6930 Big Data for Biological Applications.
The project website is at https://cis6930.wixsite.com/group3 along with the relevant literature.


## Compilation Instructions
The project is based on VARI, which can be found at : https://github.com/cosmo-team/cosmo/tree/VARI
The 3rd party packages must be installed according to the instructions on the website. The instructions there, including handling errors, are recommended. 

```sh
# Grab the code:
git clone https://github.com/theShaan/Rainbowfish-Reuse
cd Rainbowfish-Reuse/

#Follow the instructions on VARI to install the 3rd Party Dependencies. From the website -
#Note: Change "/home/Shaan/git/test/cosmo" to wherever your cosmo working tree ends up.

mkdir 3rd_party_src
mkdir -p 3rd_party_inst/boost
cd 3rd_party_src
git clone https://github.com/refresh-bio/KMC
git clone https://github.com/cosmo-team/sdsl-lite.git
git clone https://github.com/stxxl/stxxl
git clone https://github.com/eile/tclap
wget http://sourceforge.net/projects/boost/files/boost/1.54.0/boost_1_54_0.tar.bz2
tar -xjf boost_1_54_0.tar.bz2

# Build the five dependencies
cd boost_1_54_0
./bootstrap.sh --prefix=../../3rd_party_inst/boost
./b2 install
cd ..

cd sdsl-lite/
/usr/bin/time sh install.sh /home/Shaan/git/test/cosmo/3rd_party_inst
cd ..

cd stxxl
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/home/Shaan/git/test/cosmo/3rd_party_inst -DBUILD_STATIC_LIBS=ON
make
make install
cd ../..

cd KMC
make
cd ..

cd tclap/
autoreconf -fvi
./configure --prefix=/home/Shaan/git/test/cosmo/3rd_party_inst
make
make install
cd ..

# Build Solution
cd ..
make

```

## Usage

./cosmo-build <color_file> <de_bruijn_graph_file> <no_of_colors>

The output will <color_file>.reuse

Example:

To test the program on an ecoli sample run the following after compiling :

### Colored de Bruijn graph example:
```sh
# Grab 6 E. coli assemblies:
git clone https://github.com/theShaan/e_coli6
cd e_coli6

# Use KMC2 to k-mer count the FASTA (*.fna) files
mkdir kmc_temp
ls -1 --color=no *.fna |xargs -l -i  ../3rd_party_src/KMC/bin/kmc -ci0 -fm -k32 -cs300 {} {}_kmc kmc_temp
ls -1 --color=no *.fna |xargs -l -i  ../3rd_party_src/KMC/bin sort {}_kmc {}_kmc_sorted_kmc.kmc
ls -1 --color=no *.fna |xargs -l -i echo "{}_kmc_sorted_kmc.kmc" >ecoli6_kmc2_list

# Build the succinct de Bruijn graph and permute uncompresed color matrix accordingly
cosmo-build -d ecoli6_kmc2_list

# Run the reuse algorithm, the output will be the color file appended with .reuse
cosmo-build -r ecoli6_kmc2_list.colors ecoli6_kmc2_list.dbg 6
```
