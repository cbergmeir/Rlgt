rPackageName=Rlgt
newDate=$(shell date +%Y-%m-%d)
rPackageVersion=0.1-1

fixPermissions:
#sed running on Windows screws up permissions
ifeq ($(OS),Windows_NT)
	icacls $(rPackageName)/DESCRIPTION /reset
	icacls $(rPackageName)/NAMESPACE /reset
endif	

check: roxy
	R CMD build $(rPackageName)
	R CMD check $(rPackageName)_$(rPackageVersion).tar.gz

install: roxy
	R CMD build $(rPackageName)
	R CMD INSTALL $(rPackageName)_$(rPackageVersion).tar.gz

roxy:
	rm -f ./$(rPackageName)/man/*.Rd
	echo "library(roxygen2)"> tmp_roxy.R
	echo "path <- \"./$(rPackageName)\"" >> tmp_roxy.R
	echo "roxygenize(package.dir=path)" >> tmp_roxy.R
	R CMD BATCH tmp_roxy.R
	cd ./$(rPackageName) && sed -i 's/\(Date: \).*/Date: '"$(newDate)"'/' DESCRIPTION
	cd ./$(rPackageName) && sed -i -e 's/\".registration=TRUE\"/.registration=TRUE/' NAMESPACE
	
clean:
	rm -rf $(rPackageName)_$(rPackageVersion).tar.gz $(rPackageName).Rcheck
	rm -rf ./$(rPackageName)/man/*.Rd
	rm -rf ./tmp_roxy.R ./tmp_roxy.Rout
