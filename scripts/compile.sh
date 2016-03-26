#!/bin/bash
set -e

branch=$1
version=$2

echo "Building sailfish [branch = ${branch}]. Tagging version as ${version}"

# Activate Holy Build Box environment.
source /hbb_exe/activate

set -x

# Install things we need
yum install -y --quiet wget
wget http://download.fedoraproject.org/pub/epel/5/x86_64/epel-release-5-4.noarch.rpm
rpm -i --quiet epel-release-5-4.noarch.rpm
#yum install -y --quiet git
yum install -y --quiet unzip
yum install -y --quiet bzip2-devel.x86_64

curl -k -L https://github.com/kingsfordgroup/sailfish/archive/${branch}.zip -o ${branch}.zip
unzip ${branch}.zip
mv sailfish-${branch} sailfish
cd sailfish
mkdir build
cd build
cmake -DFETCH_BOOST=TRUE ..
make
make install
make test
cd ../scripts
bash make-release.sh -v ${version} -n CentOS5
cd ../RELEASES
cp *.tar.gz /io/
