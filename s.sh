for((i = 1; i <= 100 ; ++i)); do
    echo $i
    ./gen $i > in
    # ./A < int > out1
    # ./B < int > out2
    # diff -w out1 out2 || break
    # diff -w <(./A < int) <(./B < int) || break
done