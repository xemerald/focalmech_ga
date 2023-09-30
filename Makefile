#
#
#
all: libs echo_msg
	@(cd ./src; make;);
#
#
libs: echo_msg_libraries
	@(cd ./src/libsrc; make;);
#
#
echo_msg:
	@echo "---------------------------------------";
	@echo "- Making main program of focalmech_ga -";
	@echo "---------------------------------------";
echo_msg_libraries:
	@echo "----------------------------------";
	@echo "-        Making libraries        -";
	@echo "----------------------------------";

# Clean-up rules
clean:
	@(cd ./src; make clean;);
	@(cd ./src/libsrc; make clean; make clean_lib;);

clean_bin:
	@(cd ./src; make clean_bin;);
