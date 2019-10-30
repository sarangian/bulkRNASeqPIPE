#!/bin/bash
#set -e
#set -x

THIS_DIR=$(DIRNAME=$(dirname "$0"); cd "$DIRNAME"; pwd)
THIS_FILE=$(basename "$0")
THIS_PATH="$THIS_DIR/$THIS_FILE"
PREFIX=$HOME/RNASeqWF/


if [ $(which python3 2>/dev/null) ];
        then
        if [ "$(python3 --version | grep -c "3.5.")" == 1 ] || [ "$(python3 --version | grep -c "3.6.")" == 1 ]  || [ "$(python3 --version | grep -c "3.7.")" == 1 ];
                then
                echo ""
                echo -e "\e[1;36m                            Checking DISK USE ....      \e[0m"

                if [ $(which pip3 2>/dev/null) ];
                then 
                  if [ "$(pip3 --version | grep -c "19.2")" == 1 ] || [ "$(pip3 --version | grep -c "19.3")" == 1 ];
                        then
                           echo ""
                  else 
                    echo -e "\e[1;36m Installer will update $(pip3 --version) to latest \e[0m";
                  fi
                fi
        fi
fi

if [ "$(python3 --version | grep -c "3.4.")" == 1 ];
     then
     echo -e "\e[1;36m $(python3 --version) ...installed. 
     Installation of RNASeq workflow requires atleat Python 3.5 as the default version of python,
     along with pip3 minimal version 19.1.
     The installer will auto exit\e[0m";
     sleep 5s;
     exit 0
fi

sleep 2s;

echo -e "\e[1;34m__________________________________DISK USE SUMMARY______________________________\e[0m"
total=$(df --total | tail -n 1 | awk '{print $2}')
used=$(df --total | tail -n 1 | awk '{print $3}')
available=$(df --total | tail -n 1 | awk '{print $4}')

echo -e "Total Disk Space:\t $(( ${total} / 1024 /1024 )) GB"
echo -e "Used  Disk Space:\t $(( ${used} / 1024 /1024 )) GB"
echo -e "Available Disk Space:\t  $(( ${available} / 1024 /1024 )) GB"
echo -e "\e[1;34m__________________________________CPU SUMMARY____________________________________\e[0m"
lscpu | egrep 'CPU\(s\)|Thread|Model name'


echo -e "\e[1;34m__________________________________Memory USE SUMMARY_____________________________\e[0m"
totalM=$(free -m | grep "Mem:"| awk '{print $2}')
usedM=$(free -m | grep "Mem:"| awk '{print $3}')
availableM=$(free -m | grep "Mem:"| awk '{print $4}')
echo -e "Total Memory:\t $(( ${totalM} )) MB OR $(( ${totalM} / 1024 )) GB"
echo -e "Used  Memory:\t $(( ${usedM} )) MB OR $(( ${usedM} / 1024 )) GB"
echo -e "Available Memory:\t  $(( ${availableM} )) MB OR $(( ${availableM} / 1024 )) GB"

echo -e "\e[1;34m_________________________________________________________________________________\e[0m"
echo -e "Hello "$USER""
    printf "Installation RNASeq-Pipe requires at least 5gb free disk space.\\nIf you do not have sufficient disc space, Press CTRL-C to abort the installation.\\nRNASeq Workflow will now be installed into this location:\\n"
    printf "%s\\n" "$PREFIX"
    printf "\\n"
    printf "  - Press ENTER to confirm the location\\n"
    printf "  - Press CTRL-C to abort the installation\\n"
    printf "  - Or specify a different location below\\n"
    printf "\\n"
    printf "[%s] >>> " "$PREFIX"
    read -r user_prefix
    if [ "$user_prefix" != "" ]; then
        case "$user_prefix" in
            *\ * )
                printf "ERROR: Cannot install into directories with spaces\\n" >&2
                exit 1
                ;;
            *)
                eval PREFIX="$user_prefix"
                ;;
        esac
    fi

case "$PREFIX" in
    *\ * )
        printf "ERROR: Cannot install into directories with spaces\\n" >&2
        exit 1
        ;;
esac


if ! mkdir -p "$PREFIX"; then
    printf "ERROR: Could not create directory: '%s'\\n" "$PREFIX" >&2
    exit 1
fi

PREFIX=$(cd "$PREFIX"; pwd)
export PREFIX

printf "PREFIX=%s\\n" "$PREFIX"


if [ ! -d $PREFIX ]; then
  mkdir -p $PREFIX
fi
cd $PREFIX
start_dir=$(pwd)

# Make an install location
if [ ! -d 'tools' ]; then
  mkdir -p tools
fi
cd tools
build_dir=$(pwd)

######

#QC
FASTQC_VERSION=0.11.8
QUAST_VERSION=5.0.2
BBMAP_VERSION=38.61b
NCBI_BLAST_VERSION=2.9.0

#Assembly
SPADES_VERSION=3.13.0
ROCKHOPPER2_VERSION=2_0_3
CORSET_VERSION=1.09
LACE_VERSION=1.13

PROKKA_VERSION=1.14.0
PRODIGAL_VERSION=2.6.2
PPLACER_VERSION=1.1.alpha17
HMMER_VERSION=3.2


#TRANSCRIPT DE
STAR_VERSION=2.7.1a
SALMON_VERSION=0.14.1
KALLISTO_VERSION=0.46.0
HISAT2_VERSION=2.1.0
SAMTOOLS_VERSION=1.9
SUBREAD_VERSION=1.6.5

#R and PANDOC
R_VERSION=3.6.1
PANDOC_VERSION=2.7.3

#ASSEMBLY_tools_download_url
PANDOC_DOWNLOAD_URL="https://github.com/jgm/pandoc/releases/download/${PANDOC_VERSION}/pandoc-${PANDOC_VERSION}-1-amd64.deb"
FASTQC_DOWNLOAD_URL="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip"
QUAST_DOWNLOAD_URL="https://nchc.dl.sourceforge.net/project/quast/quast-${QUAST_VERSION}.tar.gz"
BBMAP_DOWNLOAD_URL="https://jaist.dl.sourceforge.net/project/bbmap/BBMap_${BBMAP_VERSION}.tar.gz"
RCORRECTOR_DOWNLOAD_URL="https://github.com/mourisl/Rcorrector.git"
RNASEQDEA_DOWNLOAD_URL="https://github.com/sarangian/RNASeqDEA.git"
CORSET_DOWNLOAD_URL="https://github.com/Oshlack/Corset/releases/download/version-${CORSET_VERSION}/corset-${CORSET_VERSION}-linux64.tar.gz"
LACE_DOWNLOAD_URL="https://github.com/Oshlack/Lace/releases/download/v1.13/Lace-${LACE_VERSION}.tar.gz"
SUBREAD_DOWNLOAD_URL="https://jaist.dl.sourceforge.net/project/subread/subread-${SUBREAD_VERSION}/subread-${SUBREAD_VERSION}-Linux-x86_64.tar.gz"
ROCKHOPPER2_DOWNLOAD_URL="https://cs.wellesley.edu/~btjaden/Rockhopper/download/version_${ROCKHOPPER2_VERSION}/Rockhopper.jar"
SPADES_DOWNLOAD_URL="http://cab.spbu.ru/files/release${SPADES_VERSION}/SPAdes-${SPADES_VERSION}-Linux.tar.gz"
#MAPPING_tools_download_url
STAR_DOWNLOAD_URL="https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz"
HISAT2_DOWNLOAD_URL="http://ccb.jhu.edu/software/hisat2/dl/hisat2-${HISAT2_VERSION}-Linux_x86_64.zip"
SALMON_DOWNLOAD_URL="https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/salmon-${SALMON_VERSION}_linux_x86_64.tar.gz"
KALLISTO_DOWNLOAD_URL="https://github.com/pachterlab/kallisto/releases/download/v${KALLISTO_VERSION}/kallisto_linux-v${KALLISTO_VERSION}.tar.gz"
DART_DOWNLOAD_URL="https://github.com/hsinnan75/Dart.git"
SAMTOOLS_DOWNLOAD_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2"

PROKKA_DOWNLOAD_URL="https://github.com/tseemann/prokka.git"
PRODIGAL_DOWNLOAD_URL="https://github.com/hyattpd/Prodigal/releases/download/v${PRODIGAL_VERSION}/prodigal.linux"
PPLACER_DOWNLOAD_URL="https://github.com/matsen/pplacer/releases/download/v${PPLACER_VERSION}/pplacer-Linux-v${PPLACER_VERSION}.zip"
#HMMER_DOWNLOAD_URL="http://eddylab.org/software/hmmer/hmmer-${HMMER_VERSION}.tar.gz"
NCBI_BLAST_DOWNLOAD_URL="ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${NCBI_BLAST_VERSION}/ncbi-blast-${NCBI_BLAST_VERSION}+-x64-linux.tar.gz"
#transcript_DE_tools_download_url


#Install R, DESEQ2, EDGER
#R_VERSION=3.6.1
release=$(lsb_release -c | cut -f2)
if [ "$release" == "xenial" ];
   then R_DOWNLOAD_URL="https://cloud.r-project.org/bin/linux/ubuntu/xenial-cran35/r-base-core_${R_VERSION}-3xenial_amd64.deb"
   #echo -e "\e[1;31m $R_DOWNLOAD_URL \e[0m";
fi

if [ "$release" == "trusty" ];
   then R_DOWNLOAD_URL="https://cloud.r-project.org/bin/linux/ubuntu/trusty-cran35/r-base-core_${R_VERSION}-3trusty2_amd64.deb"
   #echo -e "\e[1;31m $R_DOWNLOAD_URL \e[0m";
fi

if [ "$release" == "bionic" ];
   then R_DOWNLOAD_URL="https://cloud.r-project.org/bin/linux/ubuntu/bionic-cran35/r-base-core_${R_VERSION}-3bionic_amd64.deb"
   #echo -e "\e[1;31m $R_DOWNLOAD_URL \e[0m";
fi

if [ "$release" == "" ];
   then R_DOWNLOAD_URL="https://cloud.r-project.org/bin/linux/ubuntu/trusty-cran35/r-base-core_${R_VERSION}-3trusty2_amd64.deb"
   #echo -e "\e[1;31m $R_DOWNLOAD_URL \e[0m";
fi

if [ "$(R --version | grep -c "3.5.")" == 1 ] || [ "$(R --version | grep -c "3.6.")" == 1 ];
then  
    r_version=$(R --version | head -n 1 | cut -f3 -d " ")

    if [ "$r_version" == "3.6.1" ]; 
    then
    echo -e "\e[1;36m R-$r_version \t...installed \e[0m";
    else
    echo -e "\e[1;31m R version >=3.5.0 ...not in PATH \e[0m"
    fi
fi
#-----------------Checking dependancies, pre installed tools--------------------------
#----minal dependancies------
declare -a dpkglist=( "libbz2-dev" "libtcl8.6" "gfortran" "libtk8.6" "liblzma-dev" "zlib1g" "zlib1g-dev" "cpanminus" 
"libdatetime-perl" "libxml-simple-perl" "gcc" "g++" "zlib1g-dev" "liblapack-dev" "libblas-dev" "libxml2-dev" 
"libssl-dev" "libpcre3-dev" "libboost-dev" "libcurl4-gnutls-dev" "libmagick++-dev" "wget" "curl" "git" "default-jre" "pandoc"
"python" "python3" "python-docutils" "python-setuptools" "python-testresources" "python3-testresources"
"hmmer" "bowtie2" "fastqc")

## now loop through the above array
for package in "${dpkglist[@]}"
do
  if [ $(dpkg-query -W -f='${Status}' $package 2>/dev/null | grep -c "ok installed") -eq 1 ];
then
  echo -e "\e[1;36m $package \t...installed \e[0m";
else
  echo -e "\e[1;31m $package \t...not in PATH \e[0m";
fi
done



#----packages from git-hub, wget and apt-get 
#-----------------Checking pre installed tools--------------------------
declare -a tool_list=("spades.py" "quast.py" "Rockhopper.jar" "pplacer" "featureCounts" "prodigal" "samtools"
 "corset" "rcorrector" "Lace.py" "bbduk.sh" "STAR" "dart" "hisat2" "salmon" "kallisto")
## now loop through the above array
for package in "${tool_list[@]}"
do
  if [ $(which $package 2>/dev/null) ];
     then
     echo -e "\e[1;36m $package \t...installed \e[0m";
     else
     echo -e "\e[1;31m $package \t...not in PATH \e[0m";
  fi
done

#.................................BioPrl...................................
if [ "$(perl -e 'use Bio::Perl;' 2>&1 > /dev/null | grep -c "you may need to install the Bio::Perl module")" == 1 ];
    then 
    echo -e "\e[1;31m Bio::Perl \t...not in PATH \e[0m";
else
    echo -e "\e[1;36m Bio::Perl \t...installed \e[0m"
fi
echo -e "\e[1;36m Hi $USER, I will update your system and install the required tools. I will take significant time depending your internet connection speed. \e[0m";
sleep 10s;
# DOWNLOAD EXTERNAL TOOLS THROUGH WGET
download () {
  url=$1
  download_location=$2
  if [ -e $download_location ]; then
    echo -e "\e[1;36m Skipping download of $url, $download_location already exists \e[0m"
  else
    echo -e "\e[1;31m Downloading $url to $build_dir/$download_location \e[0m"
    wget $url -O $download_location
  fi
}

# DOWNLOAD EXTERNAL TOOLS FROM GITHUB USING GIT CLONE 
clone () {
  url=$1
  download_location=$2
  if [ -e $download_location ]; then
    echo -e "\e[1;36m Skipping download of $url, $download_location already exists  \e[0m"
  else
    echo -e "\e[1;31m Downloading $url to $build_dir/$download_location \e[0m"
    git clone $url $download_location
  fi
}

# INSTALLING TOOLS FROM PIP

function pip_install {
  sudo -H pip3 install "$@"
  #sudo python3 -m pip install "$@"
  if [ $? -ne 0 ]; then
    echo "could not install $p - abort"
    exit 1
  fi
}

# INSTALLING TOOLS USING APT-GET

function apt_install {
  sudo apt-get -y install $1
  if [ $? -ne 0 ]; then
    echo "could not install $1 - abort"
    exit 1
  fi
}

#-----------------------------Check and install minimal dependancies-------------------------
declare -a dpkglist=("libbz2-dev" "libtcl8.6" "gfortran" "libtk8.6" "liblzma-dev" "zlib1g" "zlib1g-dev" "cpanminus" 
"libdatetime-perl" "libxml-simple-perl" "gcc" "g++" "zlib1g-dev" "libmagick++-dev" "liblapack-dev" "libblas-dev" "libxml2-dev" 
"libssl-dev" "libpcre3-dev" "libboost-dev" "libcurl4-gnutls-dev" "wget" "curl" "git" "default-jre" 
"python" "python3" "python-docutils" "python-setuptools" "python-testresources" "python3-testresources"
"hmmer" "bowtie2" "fastqc")

## now loop through the above array
for package in "${dpkglist[@]}"
do
  if [ $(dpkg-query -W -f='${Status}' $package 2>/dev/null | grep -c "ok installed") -eq 0 ];
  then
  echo -e "installing $package"
  apt_install $package
  #sudo apt-get -f install
  #sudo apt-get -y autoremove
  #sudo apt-get update --fix-missing
  fi
done
#sudo apt-get update --fix-missing
#sudo apt-get clean

###Installing pandoc
 if [ $(dpkg-query -W -f='${Status}' $package 2>/dev/null | grep -c "ok installed") -eq 1 ];
     then
        currentver="$(pandoc --version)"
        requiredver="1.19.2.1"
        if [ "$(printf '%s\n' "$requiredver" "$currentver" | sort -V | head -n1)" = "$requiredver" ]; 
	     then 
             echo "pandoc version >=1.19.1 in path"
        else
             echo "pandoc 2.7.3 will be installed"
             cd $build_dir
             download $PANDOC_DOWNLOAD_URL "PANDOC-${PANDOC_VERSION}.deb"
             sudo dpkg -i PANDOC-${PANDOC_VERSION}.deb
        fi
else
     echo "pandoc 2.7.3 will be installed"
     cd $build_dir
     download $PANDOC_DOWNLOAD_URL "PANDOC-${PANDOC_VERSION}.deb"
     sudo dpkg -i PANDOC-${PANDOC_VERSION}.deb
fi

#########Pip3 Version
if  [python3 -c "import pip3" &> /dev/null];
     then
        currentver="$(pip3 --version)"
        requiredver="19.3.1"
        if [ "$(printf '%s\n' "$requiredver" "$currentver" | sort -V | head -n1)" = "$requiredver" ]; 
	     then 
             echo "pip version >=19.3.1 in path"
        else
             echo "pip version will be updated now"
             sudo -H pip3 install --upgrade pip
        fi
else
     echo "pip 19.3.1 will be installed"
     cd $build_dir
     pip_dir="$build_dir/pip-latest"
    		if [ ! -d "$pip_dir" ]; then
        	   mkdir -p $pip_dir
    		fi
    	         
		cd $pip_dir
    		curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
    		chmod 755 get-pip.py
    		sudo python3 get-pip.py
fi

#####
# ------------------------Check and install tools from pip --------------------------------
#if [ "$(pip3 --version | grep -c "19.")" == 1 ];
#then
    #pip_version=$(pip3 --version | cut -f2 -d " ")
    #echo -e "\e[1;36m pip-v$pip_version \tinstalled \e[0m";
#else
    #echo -e "\e[1;36m I am going to install latest version of pip \e[0m";
    #sleep 2s
    #cd $build_dir
    #pip_dir="$build_dir/pip-latest"
    #if [ ! -d "$pip_dir" ]; then
        #mkdir -p $pip_dir
    #fi
    #cd $pip_dir
    #curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
    #chmod 755 get-pip.py
    #sudo python3 get-pip.py
#fi

#----packages from pip-------
declare -a piplist=("luigi" "numpy" "scipy" "matplotlib" "pandas" "optparse" "argparse")
for package in "${piplist[@]}"
    do
      if [ "$(python3 -c "import $package" 2>&1 > /dev/null | grep -c "ModuleNotFoundError")" == 0 ];
      then
          echo -e "\e[1;36m $package \t installed \e[0m";
      else
          echo -e "\e[1;31m $package \t not in path \e[0m";
      fi
   done

#Check Biopython
if [ "$(python3 -c "import Bio" 2>&1 > /dev/null | grep -c "ModuleNotFoundError")" == 0 ];
      then
          echo -e "\e[1;36m biopython \t installed \e[0m";
      else
          echo -e "\e[1;31m biopython \t not in path \e[0m";
fi
 

declare -a piplist=("luigi" "numpy" "scipy" "matplotlib" "pandas" "optparse" "argparse")
    for package in "${piplist[@]}"
        do
            if python3 -c "import $package" &> /dev/null;
            then
            echo -e "$package \t...installed";
            else
            echo -e "\e[1;36m Installing $package  \e[0m"
            pip_install $package
            fi
    done
    
#Install Biopython
if [ "$(python3 -c "import Bio" 2>&1 > /dev/null | grep -c "ModuleNotFoundError")" == 1 ];
      then
          pip_install biopython
fi

# ----------------------------Downlaod and Install tools---------------------
if  [ "$(makeblastdb -help | grep -c "2.8.")" == 1 ] || [ "$(makeblastdb -help | grep -c "2.9.")" == 1 ];
then
        mbdb_version=$(makeblastdb -help | grep "Application to create BLAST databases" | cut -f10 -d " ")
        echo -e "\e[1;31m makeblastdb $mbdb_version available in path\e[0m";
else
        echo -e "\e[1;31m makeblastdb --version >=2.8.0 required for running PROKKA. Now downloading and installing NCBI-BLAST-2.9.0 \e[0m"; 
	download $NCBI_BLAST_DOWNLOAD_URL "BLAST-${NCBI_BLAST_VERSION}.tar.gz"
        blast_dir="$build_dir/ncbi-blast-2.9.0+/bin"
        tar -xvzf BLAST-${NCBI_BLAST_VERSION}.tar.gz
        sudo cp $blast_dir/* /usr/local/bin
fi

#-------------------------QUAST----------------------------------------------
if ! [ $(which quast.py 2 > /dev/null) ];then
    cd $build_dir
    download $QUAST_DOWNLOAD_URL "QUAST-${QUAST_VERSION}.tar.gz"
    quast_dir="$build_dir/quast-5.0.2"
    tar -xvzf QUAST-${QUAST_VERSION}.tar.gz
    cd $quast_dir
    ./install.sh
fi

#-------------------------BBMAP----------------------------------------------
if ! [ $(which bbduk.sh 2 > /dev/null) ];then
    cd $build_dir
    download $BBMAP_DOWNLOAD_URL "BBMAP-${BBMAP_VERSION}.tar.gz"
    bbmap_dir="$build_dir/bbmap"
    tar -xvzf BBMAP-${BBMAP_VERSION}.tar.gz
fi

#-------------------------RCORRECTOR------------------------------------------
if ! [ $(which rcorrector 2>/dev/null) ];then
    cd $build_dir
    clone $RCORRECTOR_DOWNLOAD_URL "R-CORRECTOR"
    rcorrector_dir="$build_dir/R-CORRECTOR"
    cd $rcorrector_dir
    make
fi

#-------------------------RNASeqDEA------------------------------------------
if ! [ $(which featureCount_edgeR.r 2>/dev/null) ];then
    cd $build_dir
    clone $RNASEQDEA_DOWNLOAD_URL "RNASEQDEA"
    rnaseqdea_dir="$build_dir/RNASEQDEA"
    chmod -R 755 $rnaseqdea_dir
fi

# --------------- --------CORSET-----------------------------------------------
if ! [ $(which corset 2>/dev/null) ];then
    cd $build_dir
    download $CORSET_DOWNLOAD_URL "corset-${CORSET_VERSION}-linux64.tar.gz"
    corset_dir="$build_dir/corset-${CORSET_VERSION}-linux64" 
    tar -xvzf corset-${CORSET_VERSION}-linux64.tar.gz
fi

# --------------- --------LACE--------------------------------------------------
if ! [ $(which Lace.py 2>/dev/null) ];then
    cd $build_dir
    download $LACE_DOWNLOAD_URL "Lace-${LACE_VERSION}.tar.gz" 
    lace_dir="$build_dir/Lace-${LACE_VERSION}" 
    tar -xvzf Lace-${LACE_VERSION}.tar.gz
fi

# --------------- -------------spades -------------------------------------------
if ! [ $(which spades.py 2>/dev/null) ];then
    cd $build_dir
    download $SPADES_DOWNLOAD_URL "SPAdes-${SPADES_VERSION}-Linux.tar.gz"
    spades_dir="$build_dir/SPAdes-${SPADES_VERSION}-Linux/bin"
    tar -zxf SPAdes-${SPADES_VERSION}-Linux.tar.gz
fi

#------------------------------Rockhopper2-------------------------------------------------
if ! [ $(which Rockhopper2 2>/dev/null) ];then
    cd $build_dir
       rockhopper_dir="$build_dir/ROCKHOPPER2"
       if [ ! -d "$rockhopper_dir" ]; then
           mkdir $rockhopper_dir
       fi
       cd $rockhopper_dir
       download $ROCKHOPPER2_DOWNLOAD_URL Rockhopper.jar
       chmod 755 Rockhopper.jar
fi

# --------------- -------------featureCounts -----------------------------------------------
if ! [ $(which featureCounts 2>/dev/null) ];then
    cd $build_dir
    download $SUBREAD_DOWNLOAD_URL "subread-${SUBREAD_VERSION}-Linux-x86_64.tar.gz"
    featureCounts_dir="$build_dir/subread-${SUBREAD_VERSION}-Linux-x86_64/bin/"
    tar -xvzf subread-${SUBREAD_VERSION}-Linux-x86_64.tar.gz
fi

# --------------- -------------pplacer ----------------------------------------------
if ! [ $(which pplacer 2>/dev/null) ];then
    cd $build_dir
    download $PPLACER_DOWNLOAD_URL "pplacer-Linux-v${PPLACER_VERSION}.zip"
    pplacer_dir="$build_dir/pplacer-Linux-v${PPLACER_VERSION}"
    unzip pplacer-Linux-v${PPLACER_VERSION}.zip
fi

# --------------- -------------prodigal ---------------------------------------------
if ! [ $(which prodigal 2>/dev/null) ];then
    cd $build_dir
    prodigal_dir="$build_dir/prodigal-${PRODIGAL_VERSION}"
       if [ ! -d "$prodigal_dir" ]; then
           mkdir $prodigal_dir
       fi
       cd $prodigal_dir
       download $PRODIGAL_DOWNLOAD_URL prodigal
       chmod 755 prodigal
fi

#-------------------------STAR Aligner----------------------------------------------
if ! [ $(which STAR 2>/dev/null) ];then  
    cd $build_dir
    download $STAR_DOWNLOAD_URL "STAR-${STAR_VERSION}-Linux.tar.gz"
    star_dir="$build_dir/STAR-${STAR_VERSION}/bin/Linux_x86_64"
    tar -xvzf STAR-${STAR_VERSION}-Linux.tar.gz
fi

#-------------------------HISAT2 Aligner----------------------------------------------
if ! [ $(which hisat2 2>/dev/null) ];then
    cd $build_dir
    download $HISAT2_DOWNLOAD_URL "hisat2-${HISAT2_VERSION}-Linux_x86_64.zip"
    hisat2_dir="$build_dir/hisat2-${HISAT2_VERSION}"
    unzip hisat2-${HISAT2_VERSION}-Linux_x86_64.zip
fi

#-------------------------DART----------------------------------------------
if ! [ $(which dart 2>/dev/null) ];then
    cd $build_dir
    clone $DART_DOWNLOAD_URL "DART"
    dart_dir="$build_dir/DART"
    cd $dart_dir
    make
fi

#-------------------------SAMTOOLS----------------------------------------------
if ! [ $(which samtools 2>/dev/null) ];then
    cd $build_dir
    download $SAMTOOLS_DOWNLOAD_URL "samtools-${SAMTOOLS_VERSION}.tar.bz2"
    tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2
    samtools_dir="$build_dir/samtools-${SAMTOOLS_VERSION}"
        cd samtools-${SAMTOOLS_VERSION}
        ./configure --prefix=$samtools_dir --without-curses
        make
        make install
fi

#-------------------------salmon----------------------------------------------
if ! [ $(which salmon 2>/dev/null) ];then
    cd $build_dir
    download $SALMON_DOWNLOAD_URL "salmon-${SALMON_VERSION}_linux_x86_64.tar.gz"
    salmon_dir="$build_dir/salmon-latest_linux_x86_64/bin"
    tar -xvzf salmon-${SALMON_VERSION}_linux_x86_64.tar.gz
fi

#-------------------------kallisto----------------------------------------------
if ! [ $(which kallisto 2>/dev/null) ];then
    cd $build_dir
    download $KALLISTO_DOWNLOAD_URL "kallisto_linux-v${KALLISTO_VERSION}.tar.gz"
    tar -xvzf kallisto_linux-v${KALLISTO_VERSION}.tar.gz
    kallisto_dir="$build_dir/kallisto"
fi

#-------------------------BIOPERL----------------------------------------------
if [ "$(perl -e 'use Bio::Perl;' 2>&1 > /dev/null | grep -c "you may need to install the Bio::Perl module")" == 1 ];
    then
    sudo cpanm Bio::Perl
else
    echo -e "Bio::Perl \t...installed";
fi

if which prokka >/dev/null; then
    echo -e "prokka \t...installed";

else
    cd $build_dir
    clone $PROKKA_DOWNLOAD_URL "PROKKA"
    prokka_dir="$build_dir/PROKKA/bin"
    cd $prokka_dir
    ./prokka --setupdb
fi


#__________________Installation of R and related packages__________________________________

if [ "$(R --version | grep -c "3.5.")" == 1 ] || [ "$(R --version | grep -c "3.6.")" == 1 ];
then
  #echo -e "\e[1;36m  R version >=3.5 ...installed \e[0m";
  r_version=$(R --version | head -n 1 | cut -f3 -d " ")
    sleep 5s;
     if [ "$(R --version | grep -c "3.5.")" == 1 ] || [ "$(R --version | grep -c "3.6.")" == 1 ];
        then
            echo -e "\e[1;36m R-$r_version \tinstalled \e[0m";

  		if [ $(sudo su - -c "R -e \"library(devtools)\"" 2>&1 | grep -c "there is no package called") -eq 1 ];
			then
  			echo -e "\e[1;31m package devtools will be installed now...\e[0m";
			sleep 2s;
			sudo su - -c "R -e \"install.packages('Rcpp')\""
			sudo su - -c "R -e \"install.packages('devtools')\""
		else
  			echo -e "\e[1;36m devtools \tinstalled \e[0m";
		fi

               
	        if [ $(sudo su - -c "R -e \"library(rnaseqdea)\"" 2>&1 | grep -c "there is no package called") -eq 1 ];
		then
  		     echo -e "\e[1;31m installing rnaseqdea \e[0m";
		     rnaseqdea_dir="$build_dir/RNASEQDEA"
		     echo -e "\e[1;31m installing rnaseqdea from $rnaseqdea_dir\e[0m";
		     sudo su - -c "R -e \"devtools::install('$rnaseqdea_dir')\""
		else
  		     echo -e "\e[1;36m RNASeqDEA \tinstalled \e[0m";
		fi
     fi	

else
  	echo -e "\e[1;31m  R is not installed. I am going to install R-3.6.1 with rnaseqdea \e[0m"
            
    sleep 2s;
    cd $build_dir
    R_dir="$build_dir/Rv361"
       if [ ! -d "$R_dir" ]; then
           mkdir -p $R_dir
       fi
       cd $R_dir
       download $R_DOWNLOAD_URL r-base-core_3.6.1-3xenial_amd64.deb
       chmod 755 r-base-core_3.6.1-3xenial_amd64.deb
       sudo dpkg -i r-base-core_3.6.1-3xenial_amd64.deb 
       sudo su - -c "R -e \"install.packages('Rcpp')\""
       sudo su - -c "R -e \"install.packages('devtools')\""
       sudo su - -c "R -e \"devtools::install('$rnaseqdea_dir')\""
fi


#--------------------------SOURCE .bashrc-------------------------------------------
#-----------------------------UPDATE PATH-------------------------------------------
echo '# Added by RNASeq Workflow  Installer' >> ~/.bashrc

cd $build_dir
update_path ()
{
    new_dir=$1
    if [[ "$PATH" =~ (^|:)"${new_dir}"(:|$) ]]
    then
        return 0
    fi
    echo 'export PATH=$PATH:'$new_dir >> ~/.bashrc

}
#......fastqc.......
#if ! [ $(which fastqc >/dev/null) ];then
#    update_path ${fastqc_dir}
#fi
#......QUAST.......
if ! [ $(which quast.py >/dev/null) ];then
    update_path ${quast_dir}
fi

#......BBDUK.......
if ! [ $(which bbduk.sh >/dev/null) ];then
    update_path ${bbmap_dir}
fi

#......Corset.......
if ! [ $(which corset >/dev/null) ];then
    update_path ${corset_dir}
fi

if ! [ $(which Lace.py >/dev/null) ];then
    update_path ${lace_dir}
fi
#......RCorrector.......
if ! [ $(which rcorrector >/dev/null) ];then
    update_path ${rcorrector_dir}
fi
#......Rockhopper2.......
if ! [ $(which Rockhopper.jar >/dev/null) ];then
    update_path ${rockhopper_dir}
fi
#......Spades.......
if ! [ $(which spades.py >/dev/null) ];then
    update_path ${spades_dir}
fi
#......pplacer.......
if ! [ $(which pplacer >/dev/null) ];then
    update_path ${pplacer_dir}
fi
#.....prodigal.......
if ! [ $(which prodigal >/dev/null) ];then
    update_path ${prodigal_dir}
fi
#.....STAR.........
if ! [ $(which STAR >/dev/null) ];then
    update_path ${star_dir}
fi

#.....DART.........
if ! [ $(which dart >/dev/null) ];then
    update_path ${dart_dir}
fi

#.....salmon.........
if ! [ $(which salmon >/dev/null) ];then
    update_path ${salmon_dir}
fi

#.....kallisto.........
if ! [ $(which kallisto >/dev/null) ];then
    update_path ${kallisto_dir}
fi

#.....hisat2.........
if ! [ $(which hisat2 >/dev/null) ];then
    update_path ${hisat2_dir}
fi

#.....samtools.........
if ! [ $(which samtools >/dev/null) ];then
    update_path ${samtools_dir}
fi

#.....featureCounts.........
if ! [ $(which featureCounts >/dev/null) ];then
    update_path ${featureCounts_dir}
fi

#.....prokka.........
if ! [ $(which prokka >/dev/null) ];then
    update_path ${prokka_dir}
fi

if ! [ $(which featureCount_DESeq2.r >/dev/null) ];then
    update_path ${rnaseqdea_dir}
fi

echo 'RNASeq Workflow Installed at' ${start_dir};
echo -e "Run \e[1;36m source ~/.bashrc \e[0m to update path";
