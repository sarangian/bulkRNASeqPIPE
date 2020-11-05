#!/usr/bin/env bash
set -e
LOGFILE=INSTALL.log
if [ -f $LOGFILE ] ; then
    rm $LOGFILE
fi

##################################################################
#Functions for progressbar and spinnors
##################################################################
spinner()
{
    local pid=$1
    local delay=4
    local spinstr='/ - \\ |'
    #local spinstr='....'
    while [ "$(ps a | awk '{print $1}' | grep $pid)" ]; do
        local temp=${spinstr#?}
        #printf "\e[97;107m . \e[0m [%c]" "$spinstr"
        printf "\e[37;47m . [%c]\e[0m" "$spinstr"
        local spinstr=$temp${spinstr%"$temp"}
        sleep $delay
        printf "\b\b\b\b\b"
    done
    printf "    \b\b\b\b"
}

##################################################################################

# Name of application to install
        echo ""
echo -e "\e[1;34mChecking dependencies for RNASeqPIPE installation ...\e[0m"  
echo ""

sleep 2s;
#Check Ubuntu
if [ -f /etc/lsb-release ]; then
    declare -a dpkglist=("make" "cmake" "g++" "gcc" "cpp" "git" "zlib1g-dev"  "python")
    for package in "${dpkglist[@]}";
	    do
  		    if [ $(dpkg-query -W -f='${Status}' $package 2>/dev/null | grep -c "ok installed") -eq 1 ];
			    then
  			    echo -e "\e[1;36m $package \t...OK...\tAvailable in PATH \e[0m";  2>&1 | tee -a $LOGFILE
		    else
  			    echo -e "\e[1;31m $package \tNot Available in PATH \e[0m";  2>&1 | tee -a $LOGFILE
			    echo -e "\e[1;33m You need to install $package using: \"sudo apt-get install $package\" \e[0m";  2>&1 | tee -a $LOGFILE 
		    fi
	    done
fi

if [ -f /etc/lsb-release ]; then
    declare -a dpkglist=("git" "zlib1g-dev" "gcc" "cpp")
    for package in "${dpkglist[@]}";
	    do
  		    if ! [ $(dpkg-query -W -f='${Status}' $package 2>/dev/null | grep -c "ok installed") -eq 1 ];
			    then
  			    exit 0
		    fi
	    done
fi


if [ -f /etc/redhat-release ]; then
    declare -a dpkglist=("git" "zlib-devel" "gcc" "gcc-c++")
    for package in "${dpkglist[@]}";
	    do
  		    if [ $(rpm -qa | grep $package 2>/dev/null | grep -c $package) -ge 1 ] ;
			    then
  			    echo -e "\e[1;36m $package \t...OK....\tAvailable in PATH \e[0m";  2>&1 | tee -a $LOGFILE
		    else
  			    echo -e "\e[1;31m $package \tNot Available in PATH \e[0m";  2>&1 | tee -a $LOGFILE
			    echo -e "\e[1;33m You need to install $package using: \"sudo yum install $package\" \e[0m";  2>&1 | tee -a $LOGFILE 
	            fi
	    done
fi

if [ -f /etc/redhat-release ]; then
    declare -a dpkglist=("git" "zlib-devel" "gcc" "gcc-c++")
    for package in "${dpkglist[@]}";
	    do
  		    if ! [ $(rpm -qa | grep $package 2>/dev/null | grep -c $package) -ge 1 ] ;
			    then
  			    exit 0
	            fi
	    done
fi

#Packages through git, wget and Sourseforge
THIS_DIR=$(DIRNAME=$(dirname "$0"); cd "$DIRNAME"; pwd)
BASE_DIR=$(pwd)
chmod 755 $BASE_DIR/*.sh
chmod 755 $BASE_DIR/*.py
echo "export PATH=\$PATH:$BASE_DIR" >> ~/.bashrc
UTILITY=$BASE_DIR/utility
DEA_SCRIPTS=$BASE_DIR/utility/deaRscripts

THIS_FILE=$(basename "$0")
THIS_PATH="$THIS_DIR/$THIS_FILE"
PREFIX=$THIS_DIR/tools

sleep 2s;

echo -e "\e[1;34m__________________________________DISK USE SUMMARY______________________________\e[0m"  2>&1 | tee -a $LOGFILE
total=$(df --total | tail -n 1 | awk '{print $2}')
used=$(df --total | tail -n 1 | awk '{print $3}')
available=$(df --total | tail -n 1 | awk '{print $4}') 

echo -e "Total Disk Space:\t $(( ${total} / 1024 /1024 )) GB"  2>&1 | tee -a $LOGFILE
echo -e "Used  Disk Space:\t $(( ${used} / 1024 /1024 )) GB"  2>&1 | tee -a $LOGFILE
echo -e "Available Disk Space:\t  $(( ${available} / 1024 /1024 )) GB"  2>&1 | tee -a $LOGFILE
echo -e "\e[1;34m__________________________________CPU SUMMARY____________________________________\e[0m"  2>&1 | tee -a $LOGFILE
lscpu | egrep 'CPU\(s\)|Thread|Model name'


echo -e "\e[1;34m__________________________________Memory USE SUMMARY_____________________________\e[0m"  2>&1 | tee -a $LOGFILE
totalM=$(free -m | grep "Mem:"| awk '{print $2}')
usedM=$(free -m | grep "Mem:"| awk '{print $3}')
availableM=$(free -m | grep "Mem:"| awk '{print $4}')
echo -e "Total Memory:\t $(( ${totalM} )) MB OR $(( ${totalM} / 1024 )) GB"  2>&1 | tee -a $LOGFILE
echo -e "Used  Memory:\t $(( ${usedM} )) MB OR $(( ${usedM} / 1024 )) GB"  2>&1 | tee -a $LOGFILE
echo -e "Available Memory:\t  $(( ${availableM} )) MB OR $(( ${availableM} / 1024 )) GB"  2>&1 | tee -a $LOGFILE

echo -e "\e[1;34m_________________________________________________________________________________\e[0m"  2>&1 | tee -a $LOGFILE
echo -e "\n"
echo -e "Hello "$USER,"" 
    printf "\nInstallation RNASeq-Pipe requires at least 5gb free disk space.\\nIf you do not have sufficient disc space, Press CTRL-C to abort the installation.\\nThird Party tools for $AppName will be installed at:\t $InstallDir"    printf "\\n"  2>&1 | tee -a $LOGFILE
    printf "%s\\n" "$PREFIX"  2>&1 | tee -a $LOGFILE
    printf "\\n"  2>&1 | tee -a $LOGFILE
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
    

AppName="RNASeqPIPE"
# Set your project's install directory name here
InstallDir=$PREFIX
#EntryPoint="YourApplicationName"
EntryPoint="RNASeqPIPE"

echo "Installing $AppName"  

sleep 2s;

# Test if new directory is empty.  Exit if it's not
if [ -d $InstallDir ]; then
    if [ "$(ls -A $InstallDir)" ]; then
        echo "ERROR: Directory is not empty" >&2
        echo "If you want to install into $InstallDir, "
        echo "clear the directory first and run this script again."
        echo "Exiting..."
        echo
        exit 1
    fi
fi


echo -e `date`  "\e[1;34m\tDownloading MiniConda... Please Wait\e[0m" 

set +e
(wget "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" -O Miniconda_Install.sh >> /dev/null 2>&1) &
spinner $!
echo -e "\n"

if [ $? -ne 0 ]; then
    (wget "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" -O Miniconda_Install.sh  >> /dev/null 2>&1) &
    spinner $!
    echo -e "\n"
fi
set -e

echo -e `date`  "\e[1;34m\tInstalling MiniConda... Please Wait\e[0m" 
#echo -e "\n"
(bash Miniconda_Install.sh -b -f -p $InstallDir >> /dev/null 2>&1) & 
spinner $!
echo -e "\n"
echo -e `date`  "\e[1;34m\tMiniConda Installed..\e[0m"
echo -e "\n"
rm Miniconda_Install.sh
# Activate the new environment
PATH="$InstallDir/bin":$PATH

# Make the new python environment completely independent
# Modify the site.py file so that USER_SITE is not imported
python -s << END
import site
site_file = site.__file__.replace(".pyc", ".py");
with open(site_file) as fin:
    lines = fin.readlines();
for i,line in enumerate(lines):
    if(line.find("ENABLE_USER_SITE = None") > -1):
        user_site_line = i;
        break;
lines[user_site_line] = "ENABLE_USER_SITE = False\n"
with open(site_file,'w') as fout:
    fout.writelines(lines)
END

#Add Script Directory

#Add Prefix to env file
#echo "prefix: $PREFIX" >> $ENVFILE
cp -ar $UTILITY $PREFIX
utility_dir=$PREFIX/utility
chmod -R 755 $utility_dir
echo "export PATH=\$PATH:$utility_dir" >> ~/.bashrc

cp -ar $DEA_SCRIPTS $PREFIX
dea_scripts_dir=$PREFIX/deaRscripts
chmod -R 755 $dea_scripts_dir
echo "export PATH=\$PATH:$dea_scripts_dir" >> ~/.bashrc
# Add Entry Point to the path
if [[ $EntryPoint ]]; then
    cd $InstallDir
    mkdir Scripts
    ln -s ../bin/$EntryPoint Scripts/$EntryPoint
    #echo "$EntryPoint script installed to $(pwd)/Scripts"  
    echo "export PATH=\"$(pwd)/Scripts\":\$PATH" >> ~/.bashrc
fi

if [[ $EntryPoint ]]; then
    cd $InstallDir
    git clone "https://github.com/grimbough/rhdf5.git" >> /dev/null 2>&1 
    rhdf5_dir="$InstallDir/rhdf5"
fi

#######
InstallDir=$PREFIX
if [[ $EntryPoint ]]; then
   cd $InstallDir
   echo -e `date` "\e[1;34m\tDownloading Salmon v1.3.0 \e[0m"
   echo -e "\n"
   (wget https://github.com/COMBINE-lab/salmon/releases/download/v1.3.0/salmon-1.3.0_linux_x86_64.tar.gz >> /dev/null 2>&1) & 
   spinner $!
   echo -e "\n"
   echo -e `date`  "\e[1;34m\tSalmon v1.3.0 Downloaded..\e[0m"
   echo -e "\n"
   tar -xvzf salmon-1.3.0_linux_x86_64.tar.gz >> /dev/null 2>&1
   mv salmon-latest_linux_x86_64/ salmon_130
   echo "export PATH=\"$InstallDir/salmon_130/bin\":\$PATH" >> ~/.bashrc

   ln -s $InstallDir/salmon_130/lib/liblzma.so.0  $InstallDir/lib/liblzma.so.0
   ln -s $InstallDir/salmon_130/lib/libm.so.6  $InstallDir/lib/libm.so.6
   ln -s $InstallDir/salmon_130/lib/libtbbmalloc_proxy.so  $InstallDir/lib/libtbbmalloc_proxy.so
   ln -s $InstallDir/salmon_130/lib/libtbbmalloc_proxy.so.2  $InstallDir/lib/libtbbmalloc_proxy.so.2
   ln -s $InstallDir/salmon_130/lib/libtbbmalloc.so  $InstallDir/lib/libtbbmalloc.so
   ln -s $InstallDir/salmon_130/lib/libtbb.so.2  $InstallDir/lib/libtbb.so.2  
fi
#######

echo -e `date`  "\e[1;34m\tCreating 'base' environment and Installing\e[0m"
declare -a list=(r-base bioconductor-delayedarray python=3.7 dart psutil libgcc imagemagick libcxx zlib cmake gxx_linux-64 cxx-compiler libxml2  libcurl libopenblas libtool curl bzip2 wget luigi pandas numpy scipy biopython perl-bioperl star  bbmap  qualimap hisat2 bowtie2 segemehl subread kallisto samtools corset lace gffread r-devtools bioconductor-summarizedexperiment  bioconductor-biocparallel bioconductor-enhancedvolcano bioconductor-genomicranges r-rcppparallel bioconductor-s4vectors bioconductor-deseq2 bioconductor-rhdf5 bioconductor-rhdf5lib r-optparse bioconductor-edger bioconductor-tximport bioconductor-genomicfeatures bioconductor-regionreport bioconductor-deformats bioconductor-plyranges r-pheatmap r-colorspace r-rcolorbrewer r-optparse r-dt r-gplots r-ggplot2 r-stringr r-tidyr r-dplyr r-rcpp r-rcpparmadillo r-readr pandoc) 
    for package in "${list[@]}";
	    do
  		#echo -e "\e[1;36m \t\t\t\tinstalling $package\e[0m"
        # 
        echo -e `date` "\e[1;36m \t\t\t\tinstalling $package\e[0m"
        echo -e  "\e[1;34mPlease Wait ...\e[0m"
        conda install -y -c anaconda -c conda-forge -c bioconda -c defaults -c statiskit -c r $package

        #(conda install -y -c anaconda -c conda-forge -c bioconda -c defaults -c statiskit -c r $package >> /dev/null 2>&1) &
        #spinner $!
        echo -e "\n"
        done

echo -e `date` "\tInstalling RHDF5 Library"
(echo "devtools::install('$rhdf5_dir')" | $InstallDir/bin/R --no-save >> /dev/null 2>&1) &
spinner $!
echo -e "\n"
echo -e `date`  "\e[1;34m\tRHDF5 Library Installed..\e[0m"
echo -e "\n"

echo -e `date` "\tInstalling rnaseqdea Scripts for Expression Analysis"
echo "devtools::install('$dea_scripts_dir')" | $InstallDir/bin/R --no-save 
echo -e "\n"
echo -e `date`  "\e[1;34m\trnaseqdea Scripts Installed..\e[0m"
echo -e "\n"

ln -s $InstallDir/lib/R/modules/lapack.so  $InstallDir/lib/libRlapack.so
ln -s $InstallDir/lib/libblas.so  $InstallDir/lib/libRblas.so
######
InstallDir=$PREFIX
if [[ $EntryPoint ]]; then
   cd $InstallDir
   (wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/v2.9.0/trinityrnaseq-v2.9.0.FULL.tar.gz  >> /dev/null 2>&1) & 
   spinner $!
   tar -xvzf trinityrnaseq-v2.9.0.FULL.tar.gz
   mv trinityrnaseq-v2.9.0 trinityrnaseq
   make -C trinityrnaseq 
   make plugins -C trinityrnaseq 
   chmod +x $InstallDir/trinityrnaseq/Trinity
   echo "export PATH=\"$InstallDir/trinityrnaseq\":\$PATH" >> ~/.bashrc
   echo "export TRINITY_HOME=$InstallDir/trinityrnaseq" >> ~/.bashrc
fi
#####

source $InstallDir/etc/profile.d/conda.sh
$InstallDir/bin/conda init bash
echo
echo -e "\e[1;36m $AppName Installed Successfully !!!!!!\e[0m"
echo ""
echo -e "\e[1;31m ****Post Installation Instructions****\e[0m"
echo -e "\e[1;36m \t[1]\tRestart the terminal first.  \e[0m" 
echo -e "\e[1;36m \t[2]\tIn the new terminal, source your .bashrc file using command: source ~/.bashrc  \e[0m" 
echo -e "\e[1;36m \t[3]\tActivate conda environment using command: conda activate \e[0m" 

echo -e "\e[1;36m Have a great day --STLab Team \e[0m"
