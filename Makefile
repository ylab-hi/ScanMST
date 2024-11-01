.PHONY: help
help: ## This help message
	@echo -e "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s :)"

clean-stubs:
	rm -rf stubs

clean: clean-stubs ## Clean up
	find .  \( -type f -name "*.py[co]" -o -type d -name "__pycache__" \) -delete && echo "Removed pycs and __pycache__"
	rm -rf dist
	rm -rf build

compile-database: ## Compile database
	bear -- poetry build

stubs: clean-stubs ## Generate pybind11 stubs
	echo "Generating pybind11 stubs"
	pybind11-stubgen scannls._cppext
	rm -rf stubs

metric:
	wily report src
