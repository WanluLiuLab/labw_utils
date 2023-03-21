.PHONY: dist
dist:
	python -m build

.PHONY: doc
doc:
	$(MAKE) -C doc

.PHONY: cleandoc
cleandoc:
	$(MAKE) -C doc clean
	$(MAKE) doc

.PHONY: serve-doc
serve-doc:
	python -m http.server -d doc/_build/html 1> /dev/null &

.PHONY: test
test:
	sh -c "PYTHONPATH=$(CURDIR)/src:$(CURDIR)/test pytest ."

.PHONY: pytype
pytype:
	 pytype --config=pytype.cfg src
