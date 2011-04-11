/*
 * linkage.h
 *
 *  Created on: 11 Apr 2011
 *      Author: benedictpaten
 */

#ifndef LINKAGE_H_
#define LINKAGE_H_

#include "cactus.h"
#include "sonLib.h"

/*
 * Picks a two points along the sequence.
 * Size of gap between the two points is picked
 * with probability log_10(x - y)/log_10(w),
 * where w is the length of the sequence.
 */
void pickAPairOfPoints(MetaSequence *metaSequence, int32_t *x, int32_t *y);

/*
 * Returns non-zero iff the two segments are within blocks that (1) both have a common sequence labelled
 * with the event identified by event string, (2) are in the same order and orientation with respect to the sequence identified
 * by (1). Both segmentX and segmentY must be part of the same (meta)sequence.
 */
bool linked(Segment *segmentX, Segment *segmentY, int32_t difference, const char *eventString);

/*
 * Samples sampleNumber pairs of positions separated by distance from 1 to bucketNumber*bucketSize.
 * int32_t *correct and *samples are arrays which must be of at least bucketNumber length, and are used
 * to accumulate sample events and record the number of correct pairs. Samples are picked using
 * pickAPairOfPoints().
 */
void samplePoints(Flower *flower, MetaSequence *metaSequence,
        const char *eventString,
        int32_t sampleNumber, int32_t *correct, int32_t *samples,
        int32_t bucketNumber, double bucketSize);

/*
 * Gets all the meta sequences in the flower that are identified by the given set of event strings.
 */
stSortedSet *getHaplotypeSequences(Flower *flower, stList *eventStrings);

#endif /* LINKAGE_H_ */
