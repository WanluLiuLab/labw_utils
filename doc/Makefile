#!/usr/bin/env make

SPHINX_OPTS := -j auto

SPHINX_APIDOC_ARGS := --separate \
	--no-toc \
	--module-first \
	--maxdepth 2 \
	--force \
	--output-dir $(CURDIR)/_apidoc \
	--doc-project 'API Docs'

BUILD_DIR := $(CURDIR)/_build
SPHINX_APIDOC := PYTHONPATH="$(dir $(CURDIR))/src:${PYTHONPATH:-}" sphinx-apidoc



.PHONY: all

all: html

.PHONY: apidoc
apidoc:
	rm -rf $(CURDIR)/_apidoc
	$(SPHINX_APIDOC) $(SPHINX_APIDOC_ARGS) $(dir $(CURDIR))/src/labw_utils
	$(SPHINX_APIDOC) $(SPHINX_APIDOC_ARGS) $(dir $(CURDIR))/src/libysjs
	$(SPHINX_APIDOC) $(SPHINX_APIDOC_ARGS) $(dir $(CURDIR))/src/ysjs
	$(SPHINX_APIDOC) $(SPHINX_APIDOC_ARGS) $(dir $(CURDIR))/src/ysjsd

.PHONY: preconfigure
preconfigure: src
	python preconfigure.py


.PHONY: src
src:
	$(MAKE) -C src


html: apidoc src refs.bibtex.bib preconfigure
	sphinx-build -M html $(CURDIR) $(BUILD_DIR) $(SPHINX_OPTS)

refs.bibtex.bib: refs.bib
	biber --tool \
		--configfile=biberConf.xml \
		--output-file refs.bibtex.bib \
		--output-legacy-date \
		--output-field-replace=location:address,journaltitle:journal \
		refs.bib

clean:
	rm -rf $(BUILD_DIR)

distclean:
