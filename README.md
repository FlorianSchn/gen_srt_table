# gen_srt_table
This program generates tables for SRT division.

## Build
### using cmake
execute
    cmake .
    make
in the respective directory
### using gcc
execute
    g++ -Wall -Wextra -Werror -std=c++1z -o gen_srt_table -O2 src/main.cpp

## Use
executing `gen_srt_table -h` shows how to use the program
    usage: ./gen_srt_table [arguments]

      -h  or --help                     print this help, other arguments are discarded
      -r <radix>                        radix a.k.a. number base in [2, 36], support for higher radices
                                        might be easy to implement, feel free (hint: add more output
                                        symbols)
                                        (default: 4)
      -d <digit range>                  digit range in ]0, radix[, the quotient digits to be used are
                                        then [-digit range, digit range]
                                        (defaults to radix-1)
      -ndb <normalized divisor bound>   the upper bound of the divisor. The usual restraint on the
                                        divisor is 1 <= divisor < r, but if the divisor is known to be
                                        smaller than another number ndb, the divisor range and therefore
                                        the table size can be decreased (defaults to radix)
      -s                                for scenarios where the partial remainder is represented as the
                                        sum of two values and is calculated after the approximation
                                        step, e.g. when carry save addition is used
      -m <max iterations>               the program tests possible SRT tables in ascending order
                                        concerning their sizes, this is the maximum number of iterations
                                        (default: 10000)
      -f <filename>                     output file, print to stdout if empty

    print configuration
      -p <print options>                various print options:
          c                             print column comments
          r                             print row comments
          h,l,a,v,u                     for cells with multiple possibilities only print the [h]ighest
                                        one, the [l]owest one, [a]ll, the highest absolute [v]alue or
                                        the lowest absolute val[u]e
                                        (defaults to h)
          i                             print table info comments
          p                             pad the partial remainder and if applicable the divisor in the
                                        comments, so every partial remainder and divisor has the same
                                        length
          t                             use two's complement for partial remainder, this can only be
                                        used if comment radix is 2 and padding with 0s is activated
                                        automatically
      -pc <comment symbol>              comment symbol to be used in the output
                                        (default: //)
      -po <oor symbol>                  symbol to use for out-of-range cells
                                        (default: 0)
      -pr <comment radix>               radix to use for the comments, has to be a multiple of radix
                                        (defaults to radix)
      -ps <separator>                   symbol to use for digit separation
                                        (default: ,)
      -pi <inner separator>             symbol to use for separation of multiple digits in the same cell
                                        (default: ,)
      -pb <begin multiple>              symbol to use for beginning multiple digits
                                        (default: [)
      -pe <end multiple>                symbol to use for ending multiple digits
                                        (default: ])

    computational hints
      -hp <p fractional>                number of fractional digits of the partial remainder including
                                        partial digits e.g. with radix 2 and -hp 2 the possible
                                        fractional digits are 00, 01, 10 and 11
      -hpp <p fractional part>          number of digits that will be skipped in the last part of the
                                        partial remainder e.g. with radix 4 and -hpp 2 possible last
                                        digits are 0 and 2; has to be in range [2, radix[ or 0 for no
                                        partial digit
      -hd <d fractional>                analogue to p fractional, but for the divisor
      -hdp <d fractional part>          analogue to p fractional part, but for the divisor
