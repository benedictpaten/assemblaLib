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

SequenceInterval *sequenceInterval_construct(int32_t start, int32_t length,
        const char *sequenceName) {
    SequenceInterval *sequenceInterval = st_malloc(sizeof(SequenceInterval));
    assert(length >= 0);
    assert(start >= 0);
    sequenceInterval->start = start;
    sequenceInterval->length = length;
    assert(sequenceName != NULL);
    sequenceInterval->sequenceName = stString_copy(sequenceName);
    return sequenceInterval;
}

void sequenceInterval_destruct(SequenceInterval *sequenceInterval) {
    free(sequenceInterval->sequenceName);
    free(sequenceInterval);
}

stList *getContigPathIntervalsP(stList *contigPaths) {
    assert(stList_length(contigPaths) > 0);
    st_logDebug("Getting contig path intervals for %i contig paths\n", stList_length(contigPaths));
    stList *intervals = stList_construct3(0, (void (*)(void *))sequenceInterval_destruct);
    for (int32_t i = 0; i < stList_length(contigPaths); i++) {
        stList *contigPath = stList_get(contigPaths, i);
        Segment *_5Segment = stList_get(contigPath, 0);
        Segment *_3Segment = stList_get(contigPath, stList_length(contigPath)-1);
        Sequence *sequence = segment_getSequence(_5Segment);

        assert(sequence != NULL);
        assert(segment_getSequence(_5Segment) == segment_getSequence(_3Segment));
        assert(segment_getStrand(_5Segment));
        assert(segment_getStrand(_3Segment));
        assert(segment_getStart(_5Segment) < segment_getStart(_3Segment) + segment_getLength(_3Segment));

        SequenceInterval *sequenceInterval = sequenceInterval_construct(segment_getStart(_5Segment),
                segment_getStart(_3Segment) + segment_getLength(_3Segment),
                sequence_getHeader(sequence));
        stList_append(intervals, sequenceInterval);
        st_logDebug("Built a path interval %s %i %i\n", sequenceInterval->sequenceName,
                sequenceInterval->start, sequenceInterval->length);
    }
    return intervals;
}

stList *getContigPathIntervals(Flower *flower, const char *chosenEventString, stList *eventStrings) {
    stList *contigPaths = getContigPaths(flower, chosenEventString, eventStrings);
    stList *intervals = getContigPathIntervalsP(contigPaths);
    stList_destruct(contigPaths);
    return intervals;
}

static int cmpFn(SequenceInterval *sequenceInteval1, SequenceInterval *sequenceInterval2) {
    assert(strcmp(sequenceInteval1->sequenceName, sequenceInterval2->sequenceName) == 0);
    assert(sequenceInteval1->start != sequenceInterval2->start);
    return sequenceInteval1->start > sequenceInterval2->start ? 1 : -1;
}

stList *getScaffoldPathIntervals(Flower *flower, const char *chosenEventString, stList *referenceEventStrings,
        stList *contaminationEventStrings, CapCodeParameters *capCodeParameters) {
    st_logDebug("Getting scaffold path intervals\n");
    stList *contigPaths = getContigPaths(flower, chosenEventString, referenceEventStrings);
    stHash *scaffoldPathsHash = getScaffoldPaths(contigPaths, referenceEventStrings, contaminationEventStrings,
            capCodeParameters);
    stList *intervals = stList_construct3(0, (void (*)(void *))sequenceInterval_destruct);
    stList *scaffoldPathsList = stHash_getValues(scaffoldPathsHash);
    stSortedSet *scaffoldPathsSet = stList_getSortedSet(scaffoldPathsList, NULL);
    stList_destruct(scaffoldPathsList);
    scaffoldPathsList = stSortedSet_getList(scaffoldPathsSet);
    for(int32_t i=0; i<stList_length(scaffoldPathsList); i++) {
        stSortedSet *scaffoldPath = stList_get(scaffoldPathsList, i);
        stList *scaffoldPathList = stSortedSet_getList(scaffoldPath);
        stList *scaffoldPathIntervals = getContigPathIntervalsP(scaffoldPathList);
        stList_sort(scaffoldPathIntervals, (int (*)(const void *, const void *))cmpFn);
        SequenceInterval *_5SequenceInterval = stList_get(scaffoldPathIntervals, 0);
        SequenceInterval *_3SequenceInterval = stList_get(scaffoldPathIntervals, stList_length(scaffoldPathIntervals)-1);
        stList_append(intervals, sequenceInterval_construct(_5SequenceInterval->start, _3SequenceInterval->start +
                _3SequenceInterval->length, _5SequenceInterval->sequenceName));
        stList_destruct(scaffoldPathList);
    }
    stList_destruct(contigPaths);
    stSortedSet_destruct(scaffoldPathsSet);
    stList_destruct(scaffoldPathsList);
    stHash_destruct(scaffoldPathsHash);
    st_logDebug("Got scaffold path intervals\n");
    return intervals;
}
