.PHONY: help
help: ## This help message
	@echo -e "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s :)"

remote: ## Sync to remote
	rsync -avhP --exclude  "*egg*" --exclude "build"  --exclude "*.so"  --exclude "poetry.lock" --exclude ".*" --exclude  "__pycache__"  ./ quest:/projects/b1171/ylk4626/project/scannls

local: ## Sync to local
	rsync -avhP  --exclude  "*egg*" --exclude "build"  --exclude "*.so"  --exclude "poetry.lock" --exclude ".*" --exclude  "__pycache__" quest:/projects/b1171/ylk4626/project/scannls ./

clean: ## Clean up
	find .  \( -type f -name "*.py[co]" -o -type d -name "__pycache__" \) -delete && echo "Removed pycs and __pycache__"
	rm -rf dist
	rm -rf build
