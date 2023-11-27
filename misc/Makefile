.PHONY: build

build:
	@echo "Unpacking CIRIquant environment"	
	conda-unpack
	sed -i '1 s/^.*$$/#!\/usr\/bin\/env python/' bin/CIRIquant
	sed -i '1 s/^.*$$/#!\/usr\/bin\/env python/' bin/CIRI_DE
	sed -i '1 s/^.*$$/#!\/usr\/bin\/env python/' bin/CIRI_DE_replicate
	sed -i '1 s/^.*$$/#!\/usr\/bin\/env python/' bin/prep_CIRIquant
	@echo "Successfully configured CIRIquant package"
