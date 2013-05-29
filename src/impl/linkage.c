/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLib.h"
#include "cactus.h"
#include "adjacencyTraversal.h"

static bool stringIsInList(const char *eventString, stList *eventStrings) {
    for (int64_t i = 0; i < stList_length(eventStrings); i++) {
        const char *eventString2 = stList_get(eventStrings, i);
        if (strcmp(eventString, eventString2) == 0) {
            return 0;
        }
    }
    return 1;
}

static void getMetaSequencesForEventsP(stSortedSet *metaSequences,
        Flower *flower, stList *eventStrings) {
    //Iterate over the sequences in the flower.
    Flower_SequenceIterator *seqIt = flower_getSequenceIterator(flower);
    Sequence *sequence;
    while ((sequence = flower_getNextSequence(seqIt)) != NULL) {
        MetaSequence *metaSequence = sequence_getMetaSequence(sequence);
        if (stringIsInList(event_getHeader(sequence_getEvent(sequence)),
                eventStrings) == 0) {
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
            getMetaSequencesForEventsP(metaSequences,
                    group_getNestedFlower(group), eventStrings);
        }
    }
    flower_destructGroupIterator(groupIt);
}

stSortedSet *getMetaSequencesForEvents(Flower *flower, stList *eventStrings) {
    /*
     * Gets the haplotype sequences in the set.
     */
    stSortedSet *metaSequences = stSortedSet_construct();
    getMetaSequencesForEventsP(metaSequences, flower, eventStrings);
    return metaSequences;
}

static void getOrderedSegmentsP(Flower *flower,
        stSortedSet *segments) {
    Flower_SegmentIterator *segmentIt = flower_getSegmentIterator(flower);
    Segment *segment;
    while ((segment = flower_getNextSegment(segmentIt)) != NULL) {
            if (!segment_getStrand(segment)) {
                segment = segment_getReverse(segment);
            }
            assert(stSortedSet_search(segments, segment) == NULL);
            stSortedSet_insert(segments, segment);
    }
    flower_destructSegmentIterator(segmentIt);
    //Recurse over the flowers
    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIt)) != NULL) {
        if (group_getNestedFlower(group) != NULL) {
            getOrderedSegmentsP(group_getNestedFlower(group),
                    segments);
        }
    }
    flower_destructGroupIterator(groupIt);
}

static int64_t segmentCompareFn_coordinate;
static Name segmentCompareFn_metaSequence;
static int segmentCompareFn(const void *segment1, const void *segment2) {
    Name name1 = segment1 == &segmentCompareFn_coordinate ? segmentCompareFn_metaSequence : metaSequence_getName(sequence_getMetaSequence(segment_getSequence((Segment *)segment1)));

    Name name2 = segment2 == &segmentCompareFn_coordinate ? segmentCompareFn_metaSequence : metaSequence_getName(sequence_getMetaSequence(segment_getSequence((Segment *)segment2)));
    int i = cactusMisc_nameCompare(name1, name2);
    if(i == 0) {
        int64_t
                x =
                        segment1 == &segmentCompareFn_coordinate ? segmentCompareFn_coordinate
                                : (int64_t) segment_getStart((void *) segment1);
        int64_t
                y =
                        segment2 == &segmentCompareFn_coordinate ? segmentCompareFn_coordinate
                                : (int64_t) segment_getStart((void *) segment2);
        assert(
                segment1 == &segmentCompareFn_coordinate || segment_getStrand(
                        (void *) segment1));
        assert(
                segment2 == &segmentCompareFn_coordinate || segment_getStrand(
                        (void *) segment2));
        return x > y ? 1 : (x < y ? -1 : 0); //i > 0 ? 1 : i < 0 ? -1 : 0; //This was because of an overflow
    }
    return i;
}

stSortedSet *getOrderedSegments(Flower *flower) {
    stSortedSet *segments = stSortedSet_construct3(segmentCompareFn, NULL);
    getOrderedSegmentsP(flower, segments);
    return segments;
}

static Segment *getSegment(stSortedSet *sortedSegments, int64_t x, MetaSequence *metaSequence) {
    segmentCompareFn_coordinate = x;
    segmentCompareFn_metaSequence = metaSequence_getName(metaSequence);
    Segment *segment = stSortedSet_searchLessThanOrEqual(sortedSegments,
            &segmentCompareFn_coordinate);
    assert((void *) segment != &segmentCompareFn_coordinate);


    if (segment != NULL) {
        if(sequence_getMetaSequence(segment_getSequence(segment)) == metaSequence) {
            assert(segment_getStart(segment) <= x);
            if (x < segment_getStart(segment) + segment_getLength(segment)) {
                return segment;
            }
        }
    }
    return NULL;
}

void pickAPairOfPointsP(MetaSequence *metaSequence, int64_t *x, int64_t *y, double proportionOfSequence) {
    assert(metaSequence_getLength(metaSequence) > 20);
    assert(proportionOfSequence > 0);
    assert(proportionOfSequence <= 1.0);
    double interval = log10(metaSequence_getLength(metaSequence) * proportionOfSequence - 10);
    double j = RANDOM();
    double i = interval * j;
    int64_t size = (int64_t) pow(10.0, i) + 1;
    assert(size >= 1);
    assert(size < metaSequence_getLength(metaSequence));
    *x = metaSequence_getStart(metaSequence) + RANDOM()
            * (metaSequence_getLength(metaSequence) - size - 5);
    *y = *x + size;
    assert(*x >= 0);
    assert(*x < *y);
    assert(
            *y - metaSequence_getStart(metaSequence) < metaSequence_getLength(
                    metaSequence));
}

void pickAPairOfPoints(MetaSequence *metaSequence, int64_t *x, int64_t *y) {
    pickAPairOfPointsP(metaSequence, x, y, 1.0);
}

bool linked(Segment *segmentX, Segment *segmentY, int64_t difference,
        const char *eventString, bool *aligned) {
    assert(segment_getStrand(segmentX));
    assert(segment_getStrand(segmentY));
    *aligned = 0;
    if (segment_getStart(segmentX) < segment_getStart(segmentY)) {
        Block *blockX = segment_getBlock(segmentX);
        Block *blockY = segment_getBlock(segmentY);
        Block_InstanceIterator *instanceItX = block_getInstanceIterator(blockX);
        Segment *segmentX2;
        while ((segmentX2 = block_getNext(instanceItX)) != NULL) {
            if (strcmp(event_getHeader(segment_getEvent(segmentX2)),
                    eventString) == 0) {
                Block_InstanceIterator *instanceItY =
                        block_getInstanceIterator(blockY);
                Segment *segmentY2;
                while ((segmentY2 = block_getNext(instanceItY)) != NULL) {
                    if (strcmp(event_getHeader(segment_getEvent(segmentY2)),
                            eventString) == 0) {
                        *aligned = 1;
                        if (sequence_getMetaSequence(
                                segment_getSequence(segmentX2))
                                == sequence_getMetaSequence(
                                        segment_getSequence(segmentY2))) { //Have the same assembly sequence
                            //Now check if the two segments are connected by a path of adjacency from the 3' end of segmentX to the 5' end of segmentY.
                            int64_t separationDistance;
                            if (capsAreAdjacent(segment_get3Cap(segmentX2),
                                    segment_get5Cap(segmentY2),
                                    &separationDistance)) {
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
    } else {
        assert(segmentX == segmentY);
        if(hasCapInEvent(block_get5End(segment_getBlock(segmentX)),
                eventString)) {
            *aligned = 1;
            return 1;
        }
    }
    return 0;
}

static bool duplicated(Segment *segment) {
    Sequence *sequence = segment_getSequence(segment);
    assert(sequence != NULL);
    MetaSequence *metaSequence = sequence_getMetaSequence(sequence);
    Block *block = segment_getBlock(segment);
    Block_InstanceIterator *it = block_getInstanceIterator(block);
    Segment *segment2;
    while((segment2 = block_getNext(it)) != NULL) {
       if(segment != segment2) {
           assert(segment != segment_getReverse(segment2));
           Sequence *sequence2 = segment_getSequence(segment2);
           if(sequence2 != NULL && sequence_getMetaSequence(sequence2) == metaSequence) {
               block_destructInstanceIterator(it);
               return 1;
           }
       }
    }
    block_destructInstanceIterator(it);
    return 0;
}

void samplePoints(Flower *flower, MetaSequence *metaSequence,
        const char *eventString, int64_t sampleNumber, int64_t *correct, int64_t *aligned,
        int64_t *samples, int64_t bucketNumber, double bucketSize, stSortedSet *sortedSegments,
        bool duplication, double proportionOfSequence) {
    if(metaSequence_getLength(metaSequence) <= 1) {
        return;
    }
    for (int64_t i = 0; i < sampleNumber; i++) {
        int64_t x, y;
        pickAPairOfPointsP(metaSequence, &x, &y, proportionOfSequence);
        int64_t diff = y - x;
        assert(diff >= 1);
        int64_t bucket = log10(diff) * bucketSize;
        //st_uglyf("I have %" PRIi64 " %f %" PRIi64 " %" PRIi64 "\n", (int64_t)diff, bucketSize, bucket, bucketNumber);
        assert(bucket < bucketNumber);
        assert(bucket >= 0);
        samples[bucket]++;
        Segment *segmentX = getSegment(sortedSegments, x, metaSequence);
        if (segmentX != NULL  && (duplication || !duplicated(segmentX))) {
            Segment *segmentY = getSegment(sortedSegments, y, metaSequence);
            if (segmentY != NULL && (duplication || !duplicated(segmentX))) {
                bool b;
                if(linked(segmentX, segmentY, diff, eventString, &b)) {
                    correct[bucket]++;
                }
                if(b) {
                    aligned[bucket]++;
                }
            }
        }
    }
}

void samplePointsWithOtherReference(Flower *flower, MetaSequence *metaSequence,
        const char *eventString, const char *otherEventString, int64_t sampleNumber, int64_t *correct, int64_t *aligned,
        int64_t *samples, int64_t bucketNumber, double bucketSize, stSortedSet *sortedSegments,
        bool duplication, double proportionOfSequence) {
    if(metaSequence_getLength(metaSequence) <= 1) {
        return;
    }
    for (int64_t i = 0; i < sampleNumber; i++) {
        int64_t x, y;
        pickAPairOfPointsP(metaSequence, &x, &y, proportionOfSequence);
        int64_t diff = y - x;
        assert(diff >= 1);
        int64_t bucket = log10(diff) * bucketSize;
        //st_uglyf("I have %" PRIi64 " %f %" PRIi64 " %" PRIi64 "\n", (int64_t)diff, bucketSize, bucket, bucketNumber);
        assert(bucket < bucketNumber);
        assert(bucket >= 0);
        samples[bucket]++;
        Segment *segmentX = getSegment(sortedSegments, x, metaSequence);
        if (segmentX != NULL  && (duplication || !duplicated(segmentX))) {
            Segment *segmentY = getSegment(sortedSegments, y, metaSequence);
            if (segmentY != NULL && (duplication || !duplicated(segmentX))) {
                bool b;
                linked(segmentX, segmentY, diff, otherEventString, &b);
                if(b) {
                    if(linked(segmentX, segmentY, diff, eventString, &b)) {
                        correct[bucket]++;
                    }
                    if(b) {
                        aligned[bucket]++;
                    }
                }
            }
        }
    }
}

