.PHONY: dist
dist: requirements_all.txt
	python -m build

.PHONY: requirements_all.txt # Fast so can always be built.
requirements_all.txt:
	cat requirements.txt \
		requirements_bioutils.txt \
		requirements_mlutils.txt \
		requirements_ysjs.txt \
		requirements_ysjsd.txt \
	| sort | uniq \
	> requirements_all.txt

.PHONY: doc
doc:
	$(MAKE) -C doc

.PHONY: twine
twine:
	twine upload \
		-r local-pypi \
		--skip-existing \
		--config-file twine.pypirc.ini \
		dist/*

.PHONY: html
html:
	$(MAKE) -C doc html

.PHONY: pdf
pdf:
	$(MAKE) -C doc pdf

.PHONY: cleandoc
cleandoc:
	$(MAKE) -C doc clean
	$(MAKE) doc

.PHONY: serve-doc
serve-doc:
	python -m http.server -d doc/_build/html

.PHONY: test
test:
	sh -c "PYTHONPATH=$(CURDIR)/src:$(CURDIR)/test pytest ."

.PHONY: pytype
pytype:
	 pytype --config=pytype.cfg src/labw_utils

.PHONY: sonar-scanner
sonar-scanner:
	sonar-scanner
