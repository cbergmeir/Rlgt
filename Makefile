rPackageName=Rlgt
newDate=$(date +%Y-%m-%d)
rPackageVersion=0.0-1

check: roxy
	R CMD build $(rPackageName)
	R CMD check $(rPackageName)_$(rPackageVersion).tar.gz

install: roxy
	R CMD build $(rPackageName)
	R CMD INSTALL $(rPackageName)_$(rPackageVersion).tar.gz

roxy:
	rm -f ./$(rPackageName)/man/*.Rd
	printf "library(roxygen2)\npath <- \"./$(rPackageName)/\"\nroxygenize(package.dir=path)\n" > tmp_roxy.R
	R CMD BATCH tmp_roxy.R
	cd ./$(rPackageName) && sed -i 's/\(Date: \).*/Date: '"$(newDate)"'/' DESCRIPTION
	cd ./$(rPackageName) && sed -i -e 's/\".registration=TRUE\"/.registration=TRUE/' NAMESPACE

clean:
	rm -rf $(rPackageName)_$(rPackageVersion).tar.gz $(rPackageName).Rcheck
	rm -rf ./$(rPackageName)/man/*.Rd
	rm -rf ./tmp_roxy.R ./tmp_roxy.Rout
