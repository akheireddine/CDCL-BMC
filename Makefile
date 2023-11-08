
##################################################
###                PaInleSS                    ###
##################################################
all:  build-maplecomsps build-painless
#	+ cd painless-src 						&& \
#	make release							&& \
#	mv painless ../painless


##################################################
###                MapleCOMSPS                 ###
##################################################
build-maplecomsps:
	if [ -d mapleCOMSPS/m4ri-20140914 ]; then : ; \
	else cd mapleCOMSPS && tar zxvf m4ri-20140914.tar.gz && \
	cd m4ri-20140914 && ./configure; fi
	+ $(MAKE) -C mapleCOMSPS/m4ri-20140914
	+ $(MAKE) -C mapleCOMSPS r
	mv mapleCOMSPS/build/release/bin/mapleCOMSPS maplecomsps

##################################################
###                 PaInleSS                   ###
##################################################
build-painless:
	+ $(MAKE) -C painless-src
	mv painless-src/painless painless-mcomsps
	
docs:
	rm -rf documents
	doxygen doxygen.config
	mkdir -p documents
	mv html latex documents

clean:
	##################################################
	###               MapleCOMSPS                  ###
	##################################################
	rm -rf mapleCOMSPS/m4ri-20140914
	+ $(MAKE) -C mapleCOMSPS clean

	##################################################
	###                 PaInleSS                   ###
	##################################################
	+ $(MAKE) clean -C painless-src
	rm -f painless-mcomsps
