remote:
   rsync -avhP --exclude "build"  --exclude "*.so"  --exclude "poetry.lock" --exclude ".*" --exclude  "__pycache__"  ./ quest:/projects/b1171/ylk4626/project/scannls

local:
   rsync -avhP  --exclude "build"  --exclude "*.so"  --exclude "poetry.lock" --exclude ".*" --exclude  "__pycache__" quest:/projects/b1171/ylk4626/project/scannls ./

clean:
  find .  \( -type f -name "*.py[co]" -o -type d -name "__pycache__" \) -delete && echo "Removed pycs and __pycache__"
  rm -rf dist
  rm -rf build


