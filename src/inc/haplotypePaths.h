/*
 * maximalHaplotypePaths.h
 *
 *  Created on: 24 Feb 2011
 *      Author: benedictpaten
 */

#ifndef MAXIMALHAPLOTYPEPATHS_H_
#define MAXIMALHAPLOTYPEPATHS_H_

#include "cactus.h"
#include "sonLib.h"



/*
 * Returns a list of haplotype paths (each haplotype path is represented by a list of segments).
 */
stList *getMaximalHaplotypePaths(Flower *flower, const char *eventString, stList *eventStrings);

/*
 * Get a hash of segments to maximal haplotypes paths.
 */
stHash *buildSegmentToMaximalHaplotypePathHash(stList *maximalHaplotypePaths);

/*
 * Returns the length of a given haplotype path.
 */
int32_t haplotypePathLength(stList *haplotypePath);

/*
 * Gets a hash of maximal haplotype paths to there lengths (represented as stIntTuples).
 */
stHash *buildMaximalHaplotypeToMaximalHaplotypeLength(
        stList *maximalHaplotypePaths);

#endif /* MAXIMALHAPLOTYPEPATHS_H_ */
