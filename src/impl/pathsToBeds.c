/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactus.h"
#include "contigPaths.h"
#include "adjacencyTraversal.h"
#include "adjacencyClassification.h"
#include "scaffoldPaths.h"
#include "adjacencyClassification.h"
#include "contigPaths.h"
#include "scaffoldPaths.h"
#include "pathsToBeds.h"

SequenceInterval *sequenceInterval_construct(int32_t start, int32_t end,
        const char *sequenceName) {
    SequenceInterval *sequenceInterval = st_malloc(sizeof(SequenceInterval));
    assert(end >= start);
    assert(start >= 0);
    sequenceInterval->start = start;
    sequenceInterval->end = end;
    assert(sequenceName != NULL);
    sequenceInterval->sequenceName = stString_copy(sequenceName);
    return sequenceInterval;
}

void sequenceInterval_destruct(SequenceInterval *sequenceInterval) {
    free(sequenceInterval->sequenceName);
    free(sequenceInterval);
}

static void addInterval(Segment *_5Segment, Segment *_3Segment,
        stList *intervals) {
    /*
     * Add an interval to the list of inserval.
     */
    assert(segment_getStrand(_5Segment) == segment_getStrand(_3Segment));
    if (!segment_getStrand(_5Segment)) {
        Segment *j = _5Segment;
        _5Segment = segment_getReverse(_3Segment);
        _3Segment = segment_getReverse(j);
    }
    Sequence *sequence = segment_getSequence(_5Segment);

    assert(sequence != NULL);
    assert(
            sequence_getMetaSequence(segment_getSequence(_5Segment))
                    == sequence_getMetaSequence(segment_getSequence(_3Segment)));
    assert(segment_getStrand(_5Segment));
    assert(segment_getStrand(_3Segment));
    assert(
            segment_getStart(_5Segment) < segment_getStart(_3Segment)
                    + segment_getLength(_3Segment));

    SequenceInterval *sequenceInterval =
            sequenceInterval_construct(
                    segment_getStart(_5Segment) - sequence_getStart(sequence),
                    segment_getStart(_3Segment) + segment_getLength(_3Segment)
                            - sequence_getStart(sequence),
                    sequence_getHeader(sequence));
    stList_append(intervals, sequenceInterval);
    st_logDebug("Built a path interval %s %i %i\n",
            sequenceInterval->sequenceName, sequenceInterval->start,
            sequenceInterval->end);
}

stList *getContigPathIntervals(Flower *flower, stList *contigPaths,
        const char *chosenEventString, stList *eventStrings) {
    assert(stList_length(contigPaths) > 0);
    st_logDebug("Getting contig path intervals for %i contig paths\n",
            stList_length(contigPaths));
    stList *intervals = stList_construct3(0,
            (void(*)(void *)) sequenceInterval_destruct);
    for (int32_t i = 0; i < stList_length(contigPaths); i++) {
        stList *contigPath = stList_get(contigPaths, i);
        addInterval(stList_get(contigPath, 0),
                stList_get(contigPath, stList_length(contigPath) - 1),
                intervals);
    }
    return intervals;
}

typedef struct _segmentAndPosition {
        Segment *segment;
        int32_t position;
} SegmentAndPosition;

static SegmentAndPosition *segmentAndPosition_construct(Segment *segment, int32_t position) {
    SegmentAndPosition *segmentAndPosition = st_malloc(sizeof(SegmentAndPosition));
    segmentAndPosition->segment = segment;
    segmentAndPosition->position = position;
    return segmentAndPosition;
}

static bool isInSet(stSortedSet *items, Segment *segment, int32_t position) {
    SegmentAndPosition *segmentAndPosition = segmentAndPosition_construct(segment, position);
    bool b = stSortedSet_search(items, segmentAndPosition) != NULL;
    free(segmentAndPosition);
    return b;
}

static void addToSet(stSortedSet *items, Segment *segment, int32_t position) {
    assert(!isInSet(items, segment, position));
    SegmentAndPosition *segmentAndPosition = segmentAndPosition_construct(segment, position);
    stSortedSet_insert(items, segmentAndPosition);
}

int segmentAndPosition_cmpFn(SegmentAndPosition *a, SegmentAndPosition *b) {
    int i = a->segment > b->segment ? 1 : a->segment < b->segment ? -1 : 0;
    if(i == 0) {
        i = a->position > b->position ? 1 : a->position < b->position ? -1 : 0;
    }
    return i;
}

static int32_t getSplitContigPathIntervalsP(Segment *segment,
        stList *contigPath, stSortedSet *seen, int32_t i) {
    addToSet(seen, segment, i);
    //st_uglyf("Gooo %lli %i %i %s %i %i\n", segment_getName(segment), segment_getStart(segment), segment_getLength(segment), sequence_getHeader(segment_getSequence(segment)), i, stList_length(contigPath));
    assert(
            segment_getBlock(segment) == segment_getBlock(
                    stList_get(contigPath, i)));
    if (i + 1 < stList_length(contigPath)) {
        Segment *segment2 = getAdjacentCapsSegment(segment_get3Cap(segment));
        if (segment2 != NULL) {
            //st_uglyf("Gooo2 %lli %i %i %s\n", segment_getName(segment2), segment_getStart(segment2), segment_getLength(segment2), sequence_getHeader(segment_getSequence(segment2)));
            assert(segment_getStrand(segment2) == segment_getStrand(segment));
            assert(getAdjacentCapsSegment(segment_get5Cap(segment2)) == segment);
            assert(!isInSet(seen, segment2, i+1));
        }
        return segment2 != NULL && segment_getBlock(segment2)
                == segment_getBlock(stList_get(contigPath, i + 1)) ? getSplitContigPathIntervalsP(
                segment2, contigPath, seen, i + 1)
                : i;
    } else {
        return i;
    }
}

static bool segmentIsInEvents(Segment *segment, stList *eventStrings) {
    for (int32_t i = 0; i < stList_length(eventStrings); i++) {
        if (strcmp(event_getHeader(segment_getEvent(segment)),
                stList_get(eventStrings, i)) == 0) {
            return 1;
        }
    }
    return 0;
}

stList *getSplitContigPathIntervals(Flower *flower, stList *contigPaths,
        const char *chosenEventString, stList *eventStrings) {
    assert(stList_length(contigPaths) > 0);
    st_logDebug("Getting split contig path intervals for %i contig paths\n",
            stList_length(contigPaths));
    stList *intervals = stList_construct3(0,
            (void(*)(void *)) sequenceInterval_destruct);
    for (int32_t i = 0; i < stList_length(contigPaths); i++) {
        stList *contigPath = stList_get(contigPaths, i);
        stSortedSet *seen = stSortedSet_construct3((int (*)(const void *, const void *))segmentAndPosition_cmpFn, free);
        assert(getSplitContigPathIntervalsP(stList_get(contigPath, 0),contigPath, seen, 0) == stList_length(contigPath)-1);
        for (int32_t j = 0; j < stList_length(contigPath); j++) {
            Segment *_5Segment = stList_get(contigPath, j);
            assert(isInSet(seen, _5Segment, j));
            Block_InstanceIterator *it = block_getInstanceIterator(
                    segment_getBlock(_5Segment));
            Segment *segment2;
            while ((segment2 = block_getNext(it)) != NULL) {
                if (segmentIsInEvents(segment2, eventStrings)) {
                    if (!isInSet(seen, segment2, j)) {
                        //st_uglyf("Starting interval\n");
                        int32_t k = getSplitContigPathIntervalsP(segment2,
                                contigPath, seen, j);
                        Segment *_3Segment = stList_get(contigPath, k);
                        addInterval(_5Segment, _3Segment, intervals);
                    }
                }
            }
            block_destructInstanceIterator(it);
        }
        stSortedSet_destruct(seen);
    }
    return intervals;
}

static int cmpFn(SequenceInterval *sequenceInteval1,
        SequenceInterval *sequenceInterval2) {
    assert(
            strcmp(sequenceInteval1->sequenceName,
                    sequenceInterval2->sequenceName) == 0);
    assert(sequenceInteval1->start != sequenceInterval2->start);
    return sequenceInteval1->start > sequenceInterval2->start ? 1 : -1;
}

stList *getScaffoldPathIntervals(Flower *flower, const char *chosenEventString,
        stList *referenceEventStrings, stList *contaminationEventStrings,
        CapCodeParameters *capCodeParameters) {
    st_logDebug("Getting scaffold path intervals\n");
    stList *contigPaths = getContigPaths(flower, chosenEventString,
            referenceEventStrings);
    stHash *scaffoldPathsHash =
            getScaffoldPaths(contigPaths, referenceEventStrings,
                    contaminationEventStrings, capCodeParameters);
    stList *intervals = stList_construct3(0,
            (void(*)(void *)) sequenceInterval_destruct);
    stList *scaffoldPathsList = stHash_getValues(scaffoldPathsHash);
    stSortedSet *scaffoldPathsSet =
            stList_getSortedSet(scaffoldPathsList, NULL);
    stList_destruct(scaffoldPathsList);
    scaffoldPathsList = stSortedSet_getList(scaffoldPathsSet);
    for (int32_t i = 0; i < stList_length(scaffoldPathsList); i++) {
        stSortedSet *scaffoldPath = stList_get(scaffoldPathsList, i);
        stList *scaffoldPathList = stSortedSet_getList(scaffoldPath);
        stList *scaffoldPathIntervals = getContigPathIntervals(flower,
                scaffoldPathList, chosenEventString, referenceEventStrings);
        stList_sort(scaffoldPathIntervals,
                (int(*)(const void *, const void *)) cmpFn);
        SequenceInterval *_5SequenceInterval = stList_get(
                scaffoldPathIntervals, 0);
        SequenceInterval *_3SequenceInterval =
                stList_get(scaffoldPathIntervals,
                        stList_length(scaffoldPathIntervals) - 1);
        stList_append(
                intervals,
                sequenceInterval_construct(_5SequenceInterval->start,
                        _3SequenceInterval->start + _3SequenceInterval->end,
                        _5SequenceInterval->sequenceName));
        stList_destruct(scaffoldPathList);
    }
    stList_destruct(contigPaths);
    stSortedSet_destruct(scaffoldPathsSet);
    stList_destruct(scaffoldPathsList);
    stHash_destruct(scaffoldPathsHash);
    st_logDebug("Got scaffold path intervals\n");
    return intervals;
}
