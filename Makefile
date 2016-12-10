newDate=$(shell date +%Y-%m-%d)

check: roxy
	R CMD build Rlgt
	R CMD check Rlgt_0.0-1.tar.gz

install: roxy
	R CMD build Rlgt
	R CMD INSTALL Rlgt_0.0-1.tar.gz

roxy:
	rm -f ./Rlgt/man/*.Rd
	cd ./Rlgt && echo -e "library(roxygen2)\npath <- \"./\"\nroxygenize(package.dir=path)\n" > tmp_roxy.R
	cd ./Rlgt && R CMD BATCH tmp_roxy.R
	cd ./Rlgt && sed -i 's/\(Date: \).*/Date: '"$(newDate)"'/' DESCRIPTION

clean:
	rm -rf Rlgt_0.0-1.tar.gz Rlgt.Rcheck
	rm -rf ./Rlgt/man/*.Rd
	rm -rf ./Rlgt/tmp_roxy.R ./Rlgt/tmp_roxy.Rout
