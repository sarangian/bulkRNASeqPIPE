#!/usr/bin/env bash
set -e

# Name of application to install

echo "# Checking dependencies for bulkRNASeqPIPE installation ..."
sleep 2s;

declare -a dpkglist=("libssl-dev" "libcurl4-openssl-dev" "libxml2-dev" "libboost-all-dev" "libbz2-dev" "liblzma-dev")
for package in "${dpkglist[@]}";
	do
  		if [ $(dpkg-query -W -f='${Status}' $package 2>/dev/null | grep -c "ok installed") -eq 1 ];
			then
  			echo -e "\e[1;36m $package \t...installed \e[0m";
		else
  			echo -e "\e[1;31m install $package manually using \"sudo apt-get install package-name\" The installer will auto exit now\e[0m";
			sleep 5s;
			exit 0
		fi
	done

#Packages through git, wget and Sourseforge

THIS_DIR=$(DIRNAME=$(dirname "$0"); cd "$DIRNAME"; pwd)
THIS_FILE=$(basename "$0")
THIS_PATH="$THIS_DIR/$THIS_FILE"
PREFIX=$HOME/RNASeqPIPE/

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
    printf "Installation RNASeq-Pipe requires at least 5gb free disk space.\\nIf you do not have sufficient disc space, Press CTRL-C to abort the installation.\\n\\nRNASeq Workflow will now be installed into this location:\\n"
    printf "\\n"
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
#cd $PREFIX

AppName="bulkRNASeqPIPE"
# Set your project's install directory name here
InstallDir=$PREFIX
#EntryPoint="YourApplicationName"
EntryPoint="bulkRNASeqPIPE"

echo
echo "Installing $AppName"

echo
echo "Installing into: $InstallDir"
echo


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

# Download and install Miniconda
set +e
curl "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" -o Miniconda_Install.sh
if [ $? -ne 0 ]; then
    curl "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" -o Miniconda_Install.sh
fi
set -e

bash Miniconda_Install.sh -b -f -p $InstallDir

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

conda config --add channels defaults
conda config --add channels anaconda
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels statiskit
conda config --add channels r
#conda install -y r-base=3.6.1 star python=3.6
conda install -y -c anaconda -c conda-forge -c bioconda -c defaults -c statiskit -c r r-base=3.6.1 r-rcurl r-devtools bioconductor-rhdf5lib gxx_linux-64 git libxml2 libcurl libopenblas libboost libtool curl bzip2 wget bbmap fastqc rcorrector spades hisat2 star corset lace salmon kallisto samtools prokka bowtie2 luigi pandas numpy scipy biopython perl-bioperl python=3.6

# Cleanup
rm Miniconda_Install.sh
conda clean -iltp --yes

# Add Entry Point to the path
if [[ $EntryPoint ]]; then

    cd $InstallDir
    mkdir Scripts
    ln -s ../bin/$EntryPoint Scripts/$EntryPoint

    echo "$EntryPoint script installed to $(pwd)/Scripts"
    echo "export PATH=\"$(pwd)/Scripts\":\$PATH" >> ~/.bashrc
fi

if [[ $EntryPoint ]]; then
    bin_dir="$InstallDir/bin"
    cd "$bin_dir"
    wget "https://cs.wellesley.edu/~btjaden/Rockhopper/download/version_2_0_3/Rockhopper.jar"
    chmod 755 Rockhopper.jar
fi

if ! [ $(which dart 2>/dev/null) ];then
    if [[ $EntryPoint ]]; then
        cd $InstallDir
        git clone "https://github.com/hsinnan75/Dart.git"
        dart_dir="$InstallDir/Dart"
        cd $dart_dir
        make
        cp bwt_index $bin_dir/
        cp dart $bin_dir/
    fi
fi

if [[ $EntryPoint ]]; then
    cd $InstallDir
    git clone https://github.com/sarangian/RNASeqDEA.git RNASEQDEA
    rnaseqdea_dir="$InstallDir/RNASEQDEA"
    chmod -R 755 $rnaseqdea_dir
    echo "export PATH=\"$rnaseqdea_dir\":\$PATH" >> ~/.bashrc
fi
#echo "devtools::install('$rnaseqdea_dir')" | $InstallDir/bin/R --no-save

source $InstallDir/etc/profile.d/conda.sh
$InstallDir/bin/conda init bash
echo
echo -e "\e[1;36m $AppName Installed Successfully !!!!!!\e[0m";
echo ""
echo -e "\e[1;31m ****Post Installation Instructions****\e[0m";
echo -e "\e[1;36m \t[1]\tRestart the terminal first.  \e[0m";
echo -e "\e[1;36m \t[2]\tIn the new terminal, source your .bashrc file using command: source ~/.bashrc  \e[0m";
echo -e "\e[1;36m \t[3]\tActivate conda environment using command: conda activate \e[0m";
echo -e "\e[1;36m Have a great day --STLab Team \e[0m";
