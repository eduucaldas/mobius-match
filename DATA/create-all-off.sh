ls non-rigid-world-DATA | grep ".vert" | cut -f 1 -d '.' | xargs -L1 ./create-off.sh

