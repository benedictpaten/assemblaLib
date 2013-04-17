/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef MAXIMALHAPLOTYPEPATHS_H_
#define MAXIMALHAPLOTYPEPATHS_H_

#include "cactus.h"
#include "sonLib.h"

/*
 * Returns a list of maximal contig paths (each contig path is represented by a list of segments).
 * The chosen event string is the identifier of the event for to get the contig paths.
 * The event strings identify the events which are used to assess contiguity of the contig paths.
 */
stList *getContigPaths(Flower *flower, const char *chosenEventString, stList *eventStrings);

/*
 * Get a hash of segments to contig paths.
 */
stHash *buildSegmentToContigPathHash(stList *contigPaths);

/*
 * Returns the length of a given contig path in terms of the total number of bases in the segments it contains.
 */
int64_t contigPathLength(stList *contigPath);

/*
 * Gets a hash of contig paths to their lengths (represented as stIntTuples).
 */
stHash *buildContigPathToContigPathLengthHash(
        stList *contigPaths);

#endif /* MAXIMALHAPLOTYPEPATHS_H_ */
