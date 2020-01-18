{
    ks[$1 $2] = $3; # save the third column using the first and second as index
    k1[$1]++;       # save the first column
    k2[$2]++;       # save the second column
}
END {                                # After processing input
    for (j in k1) {                  # loop over the first column 
        printf "\t%s", j;            # and print column headers
    };
    print "";                        # newline
    for (i in k2) {                  # loop over the second 
        printf "%s", i;              # print it as row header
        for (j in k1) {              # loop over first again
            printf "\t%s", ks[j i];  # and print values
        }
        print "";                    # newline
    }
}