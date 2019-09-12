# Make DUNE project
# Usage : project.sh 
# It uses optimized version of the DUNE.

VERSION=2.6.1
DUNE_VERSION=opt
DUNE_PATH=/opt/Dune/release-${VERSION}/opt
export DUNE_CONTROL_PATH=${DUNE_PATH}
${DUNE_PATH}/bin/duneproject

