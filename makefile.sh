#!/bin/bash
set -e
dir=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
mkdir bin
mkdir bash_scripts
mv SARS-CoV-2-*.sh SnpEff_processing.sh pi-tajima.sh ./bash_scripts/
git clone https://github.com/cfarkas/shc.git
cd shc/
./autogen.sh
./configure
make
cd ..
echo ""
echo "make done. Continue with install"
# Install

./shc/src/shc -f ./bash_scripts/SARS-CoV-2-FASTA-freebayes.sh -o ./SARS-CoV-2-FASTA-freebayes
./shc/src/shc -f ./bash_scripts/SARS-CoV-2-FASTA-freebayes-nolimit.sh -o ./SARS-CoV-2-FASTA-freebayes-nolimit
./shc/src/shc -f ./bash_scripts/SARS-CoV-2-GISAID-freebayes.sh -o ./SARS-CoV-2-GISAID-freebayes
./shc/src/shc -f ./bash_scripts/SARS-CoV-2-GISAID-freebayes-nolimit.sh -o ./SARS-CoV-2-GISAID-freebayes-nolimit
./shc/src/shc -f ./bash_scripts/SARS-CoV-2-merge-variants.sh -o ./SARS-CoV-2-merge-variants
./shc/src/shc -f ./bash_scripts/SARS-CoV-2-merge-variants-nolimit.sh -o ./SARS-CoV-2-merge-variants-nolimit
./shc/src/shc -f ./bash_scripts/SARS-CoV-2-NGS-freebayes.sh -o ./SARS-CoV-2-NGS-freebayes
./shc/src/shc -f ./bash_scripts/SARS-CoV-2-NGS-freebayes-nolimit.sh -o ./SARS-CoV-2-NGS-freebayes-nolimit
./shc/src/shc -f ./bash_scripts/SARS-CoV-2-processing-fasta.sh -o ./SARS-CoV-2-processing-fasta
./shc/src/shc -f ./bash_scripts/SARS-CoV-2-processing-fasta-nolimit.sh -o ./SARS-CoV-2-processing-fasta-nolimit
./shc/src/shc -f ./bash_scripts/SnpEff_processing.sh -o ./SnpEff_processing
./shc/src/shc -f ./bash_scripts/pi-tajima.sh -o ./pi-tajima
mv SARS-CoV-2-GISAID-freebayes SARS-CoV-2-GISAID-freebayes-nolimit SARS-CoV-2-merge-variants SARS-CoV-2-merge-variants-nolimit ./bin/
mv SARS-CoV-2-NGS-freebayes SARS-CoV-2-NGS-freebayes-nolimit SARS-CoV-2-processing-fasta SARS-CoV-2-processing-fasta-nolimit ./bin/
mv pi-tajima SnpEff_processing ./bin/
echo ""
echo "All done. Binaries are located in ./bin/ folder"
#
