/*
 * substitutions.h
 *
 *  Created on: 12 Apr 2011
 *      Author: benedictpaten
 */

#ifndef SUBSTITUTIONS_H_
#define SUBSTITUTIONS_H_

/*
 * The guess is a character from the iupac alphabet of ambiguous nucleotides.
 * If the answer is represented in the guess then the score returned
 * is equal to negative log base 2 of the number of characters the guess represents divided by 4 (bit score).
 */
double bitsScoreFn(char guess, char answer);

/*
 * Returns non-zero iff the answer is in the guess.
 */
bool correctFn(char guess, char answer);

#endif /* SUBSTITUTIONS_H_ */
