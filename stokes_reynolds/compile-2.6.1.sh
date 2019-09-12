# Compile with debug flags 
# Usage : compile.sh ProjectName [dbg]

#OPTS=dune.opts
#ROOT=$(pwd)
ROOT=/tmp
VERSION=2.6.1
DUNE_PATH=/opt/Dune/release-${VERSION}/opt     
COMMAND=all

# Broj argumenata komandne linije
argcnt=$#
if [ ${argcnt} -eq "0" ] 
then
  printf "Usage: compile.sh ProjectName [dbg] [clean]\n"
  exit 1
fi

# Ime modula
modul=$1
# Bez verzije ime direktorija je ime modula
directory=${modul}
# default je naredba all
#cmmand="all"
BUILDDIR=${ROOT}/build-release

for args in "$@"
do
	case "$args" in
		dbg) DUNE_PATH=/opt/Dune/release-${VERSION}/dbg     
                     BUILDDIR=${ROOT}/build-debug  
                     echo "Debug build!"
		     ;;
	     clean) COMMAND="make clean"
		     ;;
	esac
done

# Jedan dodatni argument
#if [ ${argcnt} -eq "2" ]
#then
#    DUNE_PATH=/opt/Dune/release-${VERSION}/dbg     
#    BUILDDIR=${ROOT}/build-debug  
#    echo "Debug build!"
#fi

export DUNE_CONTROL_PATH=${DUNE_PATH}:./${directory}
# echo  "DUNE_CONTROL_PATH="${DUNE_CONTROL_PATH}
echo "${DUNE_PATH}/bin/dunecontrol --builddir=$BUILDDIR --opts=${DUNE_PATH}/dune.opts --only=${modul} ${COMMAND}"
${DUNE_PATH}/bin/dunecontrol --builddir=$BUILDDIR --opts=${DUNE_PATH}/dune.opts --only=${modul} ${COMMAND}

