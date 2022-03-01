export PATH=~/gcc-4.8.4/bin:~/cmake-3.1.3-Linux-x86_64/bin:$PATH

cd ~

# transrate-tools
git clone --recursive https://github.com/Blahah/transrate-tools.git
cd transrate-tools/bamtools
mkdir build
cd build
cmake ..
make
cd ../..
cmake .
make
## package it up
cd src
tar zcvf bam-read_v1.0.0.beta4_linux.tar.gz bam-read
cp bam-read_v1.0.0.beta4_linux.tar.gz /vagrant/

# snap
git clone https://github.com/amplab/snap.git
cd snap
git fetch
git checkout dev
make
# package it up
tar zcvf snap_v2.0.1.tar.gz snap
cp snap_v2.0.1.tar.gz /vagrant/

# salmon
git clone https://github.com/COMBINE-lab/salmon.git
cd salmon
git fetch
git checkout develop
mkdir build
cd build
cmake -DFETCH_BOOST=TRUE -DCMAKE_INSTALL_PREFIX=~/salmon/ ..
make
make install
# collect libs for packaging
cd ..
scripts/cpld.bash bin/salmon lib
rm lib/libc.so.6
rm lib/ld-linux-x86-64l.so.2
rm lib/libdl.so.2
rm lib/libstdc++.so.6
rm lib/libgcc_s.so.1
rm lib/libpthread.so.0

# package it up
tar zcvf salmon_v1.7.0.tar.gz bin/salmon lib/
cp salmon_v1.7.0.tar.gz /vagrant/

# ruby 2.2 (for libruby)


# transrate (for the c extension)
