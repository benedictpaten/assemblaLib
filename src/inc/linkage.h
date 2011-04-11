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

bool linked(Segment *segmentX, Segment *segmentY, int32_t difference);

void samplePoints(Flower *flower, MetaSequence *metaSequence,
        const char *eventString,
        int32_t sampleNumber, int32_t *correct, int32_t *samples,
        int32_t bucketNumber, double bucketSize);

stSortedSet *getHaplotypeSequences(Flower *flower, stList *eventStrings);

#endif /* LINKAGE_H_ */
