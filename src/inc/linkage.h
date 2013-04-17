/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef LINKAGE_H_
#define LINKAGE_H_

#include "cactus.h"
#include "sonLib.h"

/*
 * Gets the segments in increasing order of the sequence.
 */
stSortedSet *getOrderedSegments(Flower *flower);

/*
 * Picks a two points along the sequence.
 * Size of gap between the two points is picked
 * with probability log_10(x - y)/log_10(w),
 * where w is the length of the sequence.
 */
void pickAPairOfPoints(MetaSequence *metaSequence, int64_t *x, int64_t *y);

/*
 * Returns non-zero iff the two segments are within blocks that (1) both have a common sequence labelled
 * with the event identified by event string, (2) are in the same order and orientation with respect to the sequence identified
 * by (1). Both segmentX and segmentY must be part of the same (meta)sequence.
 */
bool linked(Segment *segmentX, Segment *segmentY, int64_t difference, const char *eventString);

/*
 * Samples sampleNumber pairs of positions separated by distance from 1 to bucketNumber*bucketSize.
 * int64_t *correct and *samples are arrays which must be of at least bucketNumber length, and are used
 * to accumulate sample events and record the number of correct pairs. The aligned array
 * records samples which are aligned, but not neccesarily correctly linked. Samples are picked using
 * pickAPairOfPoints().
 */
void samplePoints(Flower *flower, MetaSequence *metaSequence,
        const char *eventString,
        int64_t sampleNumber, int64_t *correct, int64_t *aligned, int64_t *samples,
        int64_t bucketNumber, double bucketSize, stSortedSet *sortedSegments,  bool duplication, double proportionOfSequence);

/*
 * Gets all the meta sequences in the flower that are identified by the given set of event strings.
 */
stSortedSet *getMetaSequencesForEvents(Flower *flower, stList *eventStrings);

#endif /* LINKAGE_H_ */
