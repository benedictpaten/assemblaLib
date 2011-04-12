/*
 * substitutions.c
 *
 *  Created on: 12 Apr 2011
 *      Author: benedictpaten
 */

#include <ctype.h>
#include "sonLib.h"
#include "cactus.h"

double bitsScoreFn(char guess, char answer) {
    guess = toupper(guess);
    answer = toupper(answer);
    double f = 0.415037499278844; //-log_2(3/4)
    assert(answer == 'A' || answer == 'C' || answer == 'G' || answer == 'T');
    switch (guess) {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
            return guess == answer ? 2 : 0;
        case 'W':
            return (answer == 'A' || answer == 'T') ? 1 : 0;
        case 'S':
            return (answer == 'C' || answer == 'G') ? 1 : 0;
        case 'M':
            return (answer == 'A' || answer == 'C') ? 1 : 0;
        case 'K':
            return (answer == 'G' || answer == 'T') ? 1 : 0;
        case 'R':
            return (answer == 'A' || answer == 'G') ? 1 : 0;
        case 'Y':
            return (answer == 'C' || answer == 'T') ? 1 : 0;
        case 'B':
            return (answer != 'A') ? f : 0;
        case 'D':
            return (answer != 'C') ? f : 0;
        case 'H':
            return (answer != 'G') ? f : 0;
        case 'V':
            return (answer != 'T') ? f : 0;
        case 'N':
            return 0;
        default:
            st_uglyf("I have %c %c\n", guess, answer);
            assert(0);
    }
}

bool correctFn(char guess, char answer) {
    return (bitsScoreFn(guess, answer) != 0) || toupper(guess) == 'N';
}
