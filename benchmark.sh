if [ -f "nbody" ]; then
    for threads in 1 2 4 8 16 32
    do
        echo "+-----------------------+"
        echo "Execution time using $threads threads(s), 5000 randomly generated bodies, and 1 iteration"
        time ./nbody 1 100 -b 5000 $threads
        echo ""
    done
else
    echo "Please execute 'bash make.sh' before executing this script"
fi