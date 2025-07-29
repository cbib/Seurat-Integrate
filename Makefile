MAKEFLAGS = --silent --no-print-directory

PKG_NAME != grep -i "^package" DESCRIPTION | sed -E 's/^package[:[:space:]]+//I'
PKG_VER != grep -i "^version" DESCRIPTION | sed -E 's/^version[:[:space:]]+//I'

.PHONY: rcpp document build install full-install check website

rcpp:
	Rscript -e 'Rcpp::compileAttributes(pkgdir = ".", verbose = TRUE)'

document:
	Rscript -e 'devtools::document(pkg = ".")'

build:
	Rscript -e 'devtools::build(pkg = ".", vignettes = FALSE)'

install:
	Rscript -e 'remotes::install_local(path = "../$(PKG_NAME)_$(PKG_VER).tar.gz", upgrade = "never", repos = BiocManager::repositories())'

full-install:
	${MAKE} rcpp
	${MAKE} document
	${MAKE} build
	${MAKE} install

check:
	Rscript -e 'devtools::check(pkg = ".", vignettes = FALSE)'

website:
	Rscript -e 'pkgdown::build_site(pkg = ".")'

test:
	Rscript -e 'devtools::test(pkg = ".")'
