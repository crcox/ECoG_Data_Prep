export MCC=$(MATLABDIR)/bin/mcc
export MEX=$(MATLABDIR)/bin/mex
export MFLAGS=-m -R -singleCompThread -R -nodisplay -R -nojvm
TOP := $(shell pwd)
SRCTAR=source_code.tar.gz
SRC=.
INCL= -N -p stats -I $(SRC)
.PHONEY: all clean-all clean-postbuild glmnet sdist

all: ECOG_setup_data clean-postbuild

ECOG_setup_data: $(SRC)/setup_data.m
	$(MCC) -v $(MFLAGS) $(INCL) -o $@ $<

clean-postbuild:
	-rm *.dmr
	-rm mccExcludedFiles.log
	-rm readme.txt

sdist:
	git archive --format tar HEAD | gzip -9 > $(SRCTAR)

clean-all:
	-rm ECOG_setup_data
	-rm requiredMCRProducts.txt
	-rm build-csf-openingwindow_rsa_ecog.sh.e*
	-rm build-csf-openingwindow_rsa_ecog.sh.o*
