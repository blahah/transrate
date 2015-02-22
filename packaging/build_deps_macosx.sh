brew tap homebrew/dupes
brew tap homebrew/versions
brew install cmake zlib gcc49 autoconf

cd ~


# transrate-tools
git clone --recursive https://github.com/Blahah/transrate-tools.git
## build bamtools
cd transrate-tools/bamtools
mkdir build
cd build
cmake ..
make
cd ../..
## build transrate-tools, directly specifying the static LIBZ location
cmake -DZLIB_LIBRARY=/usr/local/opt/zlib/lib/libz.a .
make
## package it up
cd src
tar zcvf bam-read_v1.0.0.beta4_osx.tar.gz bam-read
cp bam-read_v1.0.0.beta4_osx.tar.gz /vagrant/


# snap
cd ~
git clone https://github.com/Blahah/snap.git
cd snap
git fetch
git checkout dev
mv Makefile Makefile.old
# make sure snap uses gcc
sed 's/g\+\+/g\+\+-4\.9/g' Makefile.old > Makefile
make
# package it up
tar zcvf snap_v1.0dev.67.trfix1.tar.gz snap
cp snap_v1.0dev.67.trfix1.tar.gz /vagrant/


# salmon
git clone https://github.com/kingsfordgroup/sailfish.git
cd sailfish
git fetch
git checkout develop
mkdir build
cd build
cmake -DFETCH_BOOST=TRUE ..
make
make install
cd ..
# package it up
install_name_tool -add_rpath ../lib bin/salmon
tar zcvf salmon_v0.3.0.tar.gz bin/salmon lib/
cp salmon_v0.3.0.tar.gz /vagrant/
