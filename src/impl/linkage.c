/*
 * copyNumberStats.c
 *
 *
 * Calculates stats on the different copy number of columns
 *
 *  Created on: 21 Mar 2011
 *      Author: benedictpaten
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include <ctype.h>
#include <math.h>

#include "sonLib.h"
#include "cactus.h"
#include "cactusMafs.h"
#include "assemblyStructures.h"

static void getHaplotypeSequencesP(stSortedSet *metaSequences, Flower *flower, stList *eventStrings) {
    //Iterate over the sequences in the flower.
    Flower_SequenceIterator *seqIt = flower_getSequenceIterator(flower);
    Sequence *sequence;
    while ((sequence = flower_getNextSequence(seqIt)) != NULL) {
        MetaSequence *metaSequence = sequence_getMetaSequence(sequence);
        if (strcmp(event_getHeader(sequence_getEvent(sequence)), "hapA1") == 0
                || strcmp(event_getHeader(sequence_getEvent(sequence)), "hapA2")
                        == 0) {
            if (stSortedSet_search(metaSequences, metaSequence) == NULL) {
                stSortedSet_insert(metaSequences, metaSequence);
            }
        }
    }
    flower_destructSequenceIterator(seqIt);
    //Recurse over the flowers
    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIt)) != NULL) {
        if (group_getNestedFlower(group) != NULL) {
            getHaplotypeSequencesP(metaSequences, group_getNestedFlower(group), eventStrings);
        }
    }
    flower_destructGroupIterator(groupIt);
}

stSortedSet *getHaplotypeSequences(Flower *flower, stList *eventStrings) {
    /*
     * Gets the haplotype sequences in the set.
     */
    stSortedSet *metaSequences = stSortedSet_construct();
    getHaplotypeSequencesP(metaSequences, flower, eventStrings);
    return metaSequences;
}

static void getOrderedSegmentsForSequenceP(Flower *flower,
        MetaSequence *metaSequence, stSortedSet *segments) {
    Flower_SegmentIterator *segmentIt = flower_getSegmentIterator(flower);
    Segment *segment;
    while ((segment = flower_getNextSegment(segmentIt)) != NULL) {
        if (sequence_getMetaSequence(segment_getSequence(segment))
                == metaSequence) {
            if (!segment_getStrand(segment)) {
                segment = segment_getReverse(segment);
            }
            assert(stSortedSet_search(segments, segment) == NULL);
            stSortedSet_insert(segments, segment);
        }
    }
    flower_destructSegmentIterator(segmentIt);
    //Recurse over the flowers
    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIt)) != NULL) {
        if (group_getNestedFlower(group) != NULL) {
            getOrderedSegmentsForSequenceP(group_getNestedFlower(group),
                    metaSequence, segments);
        }
    }
    flower_destructGroupIterator(groupIt);
}

static int32_t segmentCompareFn_coordinate;
static int segmentCompareFn(const void *segment1, const void *segment2) {
    int32_t
            x =
                    segment1 == &segmentCompareFn_coordinate ? segmentCompareFn_coordinate
                            : (int64_t) segment_getStart((void *) segment1);
    int32_t
            y =
                    segment2 == &segmentCompareFn_coordinate ? segmentCompareFn_coordinate
                            : (int64_t) segment_getStart((void *) segment2);
    assert(segment1 == &segmentCompareFn_coordinate || segment_getStrand(
            (void *) segment1));
    assert(segment2 == &segmentCompareFn_coordinate || segment_getStrand(
            (void *) segment2));
    int64_t i = ((int64_t) x) - y;
    return i > 0 ? 1 : i < 0 ? -1 : 0; //This was because of an overflow
}

static stSortedSet *getOrderedSegmentsForSequence(Flower *flower,
        MetaSequence *metaSequence) {
    /*
     * Gets the segments in increasing order of the sequence.
     */
    stSortedSet *segments = stSortedSet_construct3(segmentCompareFn, NULL);
    getOrderedSegmentsForSequenceP(flower, metaSequence, segments);
    //Debug check
    stSortedSetIterator *it = stSortedSet_getIterator(segments);
    Segment *segment;
    int32_t i = metaSequence_getStart(metaSequence);
    while ((segment = stSortedSet_getNext(it)) != NULL) {
        assert(i == segment_getStart(segment));
        i += segment_getLength(segment);
    }
    assert(i - metaSequence_getStart(metaSequence) == metaSequence_getLength(
            metaSequence));
    stSortedSet_destructIterator(it);
    return segments;
}

static Segment *getSegment(stSortedSet *sortedSegments, int32_t x) {
    segmentCompareFn_coordinate = x;
    Segment *segment = stSortedSet_searchLessThanOrEqual(sortedSegments,
            &segmentCompareFn_coordinate);
    assert((void *) segment != &segmentCompareFn_coordinate);
    assert(segment != NULL);
    return segment;
}

void pickAPairOfPoints(MetaSequence *metaSequence, int32_t *x, int32_t *y) {
    assert(metaSequence_getLength(metaSequence) > 20);
    double interval = log10(metaSequence_getLength(metaSequence)-10);
    double j = RANDOM();
    double i = interval * j;
    int32_t size = (int32_t) pow(10.0, i) + 1;
    assert(size >= 1);
    assert(size < metaSequence_getLength(metaSequence));
    *x = metaSequence_getStart(metaSequence) + RANDOM()
            * (metaSequence_getLength(metaSequence) - size - 5);
    *y = *x + size;
    assert(*x >= 0);
    assert(*x < *y);
    assert(*y - metaSequence_getStart(metaSequence) < metaSequence_getLength(
            metaSequence));
}

bool linked(Segment *segmentX, Segment *segmentY, int32_t difference,
        const char *eventString) {
    assert(segment_getStrand(segmentX));
    assert(segment_getStrand(segmentY));
    if(segment_getStart(segmentX) < segment_getStart(segmentY)) {
        Block *blockX = segment_getBlock(segmentX);
        Block *blockY = segment_getBlock(segmentY);
        Block_InstanceIterator *instanceItX = block_getInstanceIterator(blockX);
        Segment *segmentX2;
        while ((segmentX2 = block_getNext(instanceItX)) != NULL) {
            if (strcmp(event_getHeader(segment_getEvent(segmentX2)), eventString)
                    == 0) {
                Block_InstanceIterator *instanceItY = block_getInstanceIterator(
                        blockY);
                Segment *segmentY2;
                while ((segmentY2 = block_getNext(instanceItY)) != NULL) {
                    if (strcmp(event_getHeader(segment_getEvent(segmentY2)),
                            eventString) == 0) {
                        if (sequence_getMetaSequence(segment_getSequence(segmentX2))
                                == sequence_getMetaSequence(segment_getSequence(
                                        segmentY2))) { //Have the same assembly sequence
                            //Now check if the two segments are connected by a path of adjacency from the 3' end of segmentX to the 5' end of segmentY.
                            int32_t separationDistance;
                            if (capsAreAdjacent(segment_get3Cap(segmentX2),
                                    segment_get5Cap(segmentY2), &separationDistance)) {
                                //if(difference < 10000 || (separationDistance <=  difference * 1.5 && difference <= separationDistance * 1.5)) {
                                    block_destructInstanceIterator(instanceItX);
                                    block_destructInstanceIterator(instanceItY);
                                    return 1;
                                //}
                            }
                        }
                    }
                }
                block_destructInstanceIterator(instanceItY);
            }
        }
        block_destructInstanceIterator(instanceItX);
    }
    else {
        assert(segmentX == segmentY);
        return hasCapInEvent(block_get5End(segment_getBlock(segmentX)), eventString); //isAssemblyEnd(block_get5End(segment_getBlock(segmentX)));
    }
    return 0;
}

void samplePoints(Flower *flower, MetaSequence *metaSequence,
        const char *eventString,
        int32_t sampleNumber, int32_t *correct, int32_t *samples,
        int32_t bucketNumber, double bucketSize) {
    stSortedSet *sortedSegments = getOrderedSegmentsForSequence(flower,
            metaSequence);
    assert(metaSequence_getLength(metaSequence) > 1);
    for (int32_t i = 0; i < sampleNumber; i++) {
        int32_t x, y;
        pickAPairOfPoints(metaSequence, &x, &y);
        int64_t diff = y - x;
        assert(diff >= 1);
        int32_t bucket = log10(diff) * bucketSize;
        //st_uglyf("I have %i %f %i %i\n", (int32_t)diff, bucketSize, bucket, bucketNumber);
        assert(bucket < bucketNumber);
        assert(bucket >= 0);
        samples[bucket]++;
        Segment *segmentX = getSegment(sortedSegments, x);
        Segment *segmentY = getSegment(sortedSegments, y);
        if (linked(segmentX, segmentY, diff, eventString)) {
            correct[bucket]++;
        }
    }
    stSortedSet_destruct(sortedSegments);
}


