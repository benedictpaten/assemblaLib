#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include <ctype.h>

#include "sonLib.h"
#include "cactus.h"
#include "assemblyStructures.h"
#include "haplotypePaths.h"


static void getMaximalHaplotypePathsP3(Segment *segment,
        stList *maximalHaplotypePath, stSortedSet *segmentSet, stList *eventStrings) {
    stList_append(maximalHaplotypePath, segment);
    assert(stSortedSet_search(segmentSet, segment) == NULL);
    assert(stSortedSet_search(segmentSet, segment_getReverse(segment)) == NULL);
    stSortedSet_insert(segmentSet, segment);
    Cap *_3Cap = segment_get3Cap(segment);
    if (trueAdjacency(_3Cap, eventStrings)) { //Continue on..
        Segment *otherSegment = getAdjacentCapsSegment(_3Cap);
        if (otherSegment != NULL) {
            getMaximalHaplotypePathsP3(otherSegment, maximalHaplotypePath,
                    segmentSet, eventStrings);
        }
    }
}

static void getMaximalHaplotypePathsP2(Segment *segment,
        stList *maximalHaplotypePath, stSortedSet *segmentSet, stList *eventStrings) {
    /*
     * Iterate all the way to one end of the contig then start the traversal to define the maximal
     * haplotype path.
     */
    Cap *_5Cap = segment_get5Cap(segment);
    assert(hasCapInEvents(cap_getEnd(segment_get3Cap(segment)), eventStrings)); //isHaplotypeEnd(cap_getEnd(segment_get3Cap(segment))));
    if (trueAdjacency(_5Cap, eventStrings)) { //Check that the adjacency is supported by a haplotype path
        Segment *otherSegment = getAdjacentCapsSegment(_5Cap);
        assert(segment != otherSegment);
        assert(segment_getReverse(segment) != otherSegment);
        if (otherSegment != NULL) {
            assert(stSortedSet_search(segmentSet, otherSegment) == NULL);
            assert(stSortedSet_search(segmentSet, segment_getReverse(
                    otherSegment)) == NULL);
            assert(hasCapInEvents(cap_getEnd(segment_get3Cap(otherSegment)), eventStrings)); //isHaplotypeEnd(cap_getEnd(segment_get3Cap(otherSegment))));
            getMaximalHaplotypePathsP2(otherSegment, maximalHaplotypePath,
                    segmentSet, eventStrings);
        } else { //We need to start the maximal haplotype recursion
            getMaximalHaplotypePathsP3(segment, maximalHaplotypePath,
                    segmentSet, eventStrings);
        }
    } else {
        getMaximalHaplotypePathsP3(segment, maximalHaplotypePath, segmentSet, eventStrings);
    }
}

static void getMaximalHaplotypePathsP(Flower *flower,
        stList *maximalHaplotypePaths, stSortedSet *segmentSet,
        const char *eventString,
        stList *eventStrings) {
    /*
     *  Iterate through the segments in this flower.
     */
    Flower_SegmentIterator *segmentIt = flower_getSegmentIterator(flower);
    Segment *segment;
    while ((segment = flower_getNextSegment(segmentIt)) != NULL) {
        if (stSortedSet_search(segmentSet, segment) == NULL
                && stSortedSet_search(segmentSet, segment_getReverse(segment))
                        == NULL) { //Check we haven't yet seen this segment
            if (strcmp(event_getHeader(segment_getEvent(segment)), eventString)
                    == 0) { //Check if the segment is in the assembly
                if (hasCapInEvents(cap_getEnd(segment_get5Cap(segment)), eventStrings)) { //Is a block in a haplotype segment
                    assert(hasCapInEvents(cap_getEnd(segment_get3Cap(segment)), eventStrings)); //isHaplotypeEnd(cap_getEnd(segment_get3Cap(segment))));
                    stList *maximalHaplotypePath = stList_construct();
                    stList_append(maximalHaplotypePaths, maximalHaplotypePath);
                    getMaximalHaplotypePathsP2(segment, maximalHaplotypePath,
                            segmentSet, eventStrings);
                } else {
                    assert(!hasCapInEvents(cap_getEnd(segment_get3Cap(segment)), eventStrings));//assert(!isHaplotypeEnd(cap_getEnd(segment_get3Cap(segment))));
                }
            }
        }
    }
    flower_destructSegmentIterator(segmentIt);
    /*
     * Now recurse on the contained flowers.
     */
    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIt)) != NULL) {
        if (group_getNestedFlower(group) != NULL) {
            getMaximalHaplotypePathsP(group_getNestedFlower(group),
                    maximalHaplotypePaths, segmentSet, eventString, eventStrings);
        }
    }
    flower_destructGroupIterator(groupIt);
}

static void getMaximalHaplotypePathsCheck(Flower *flower,
        stSortedSet *segmentSet, const char *eventString, stList *eventStrings) {
    /*
     * Do debug checks that the haplotypes paths are well formed.
     */
    Flower_SegmentIterator *segmentIt = flower_getSegmentIterator(flower);
    Segment *segment;
    while ((segment = flower_getNextSegment(segmentIt)) != NULL) {
        if (strcmp(event_getHeader(segment_getEvent(segment)), "assembly") == 0) {
            if (hasCapInEvents(cap_getEnd(segment_get5Cap(segment)), eventStrings)) { //isHaplotypeEnd(cap_getEnd(segment_get5Cap(segment)))) {
                assert(stSortedSet_search(segmentSet, segment) != NULL
                        || stSortedSet_search(segmentSet, segment_getReverse(
                                segment)) != NULL);
            }
        }
    }
    flower_destructSegmentIterator(segmentIt);

    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIt)) != NULL) {
        if (group_getNestedFlower(group) != NULL) {
            getMaximalHaplotypePathsCheck(group_getNestedFlower(group),
                    segmentSet, eventString, eventStrings);
        }
    }
    flower_destructGroupIterator(groupIt);
}

stList *getContigPaths(Flower *flower, const char *eventString, stList *eventStrings) {
    stList *maximalHaplotypePaths = stList_construct3(0,
            (void(*)(void *)) stList_destruct);
    stSortedSet *segmentSet = stSortedSet_construct();
    getMaximalHaplotypePathsP(flower, maximalHaplotypePaths, segmentSet, eventString, eventStrings);

    //Do some debug checks..
    st_logDebug("We have %i maximal haplotype paths\n", stList_length(
            maximalHaplotypePaths));
    getMaximalHaplotypePathsCheck(flower, segmentSet, eventString, eventStrings);
    for (int32_t i = 0; i < stList_length(maximalHaplotypePaths); i++) {
        stList *maximalHaplotypePath = stList_get(maximalHaplotypePaths, i);
        st_logDebug("We have a maximal haplotype path with length %i\n",
                stList_length(maximalHaplotypePath));
        assert(stList_length(maximalHaplotypePath) > 0);
        Segment *_5Segment = stList_get(maximalHaplotypePath, 0);
        Segment *_3Segment = stList_get(maximalHaplotypePath, stList_length(
                maximalHaplotypePath) - 1);
        if (getAdjacentCapsSegment(segment_get5Cap(_5Segment)) != NULL) {
            assert(!trueAdjacency(segment_get5Cap(_5Segment), eventStrings));
        }
        if (getAdjacentCapsSegment(segment_get3Cap(_3Segment)) != NULL) {
            assert(!trueAdjacency(segment_get3Cap(_3Segment), eventStrings));
        }
        for (int32_t j = 0; j < stList_length(maximalHaplotypePath) - 1; j++) {
            _5Segment = stList_get(maximalHaplotypePath, j);
            _3Segment = stList_get(maximalHaplotypePath, j + 1);
            assert(trueAdjacency(segment_get3Cap(_5Segment), eventStrings));
            assert(trueAdjacency(segment_get5Cap(_3Segment), eventStrings));
            assert(cap_getAdjacency(getTerminalCap(segment_get3Cap(_5Segment)))
                    == getTerminalCap(segment_get5Cap(_3Segment)));
            assert(strcmp(event_getHeader(segment_getEvent(_5Segment)),
                   eventString) == 0);
            assert(strcmp(event_getHeader(segment_getEvent(_3Segment)),
                    eventString) == 0);
            assert(hasCapInEvents(cap_getEnd(segment_get5Cap(_5Segment)), eventStrings)); //isHaplotypeEnd(cap_getEnd(segment_get5Cap(_5Segment))));
            assert(hasCapInEvents(cap_getEnd(segment_get5Cap(_3Segment)), eventStrings)); //isHaplotypeEnd(cap_getEnd(segment_get5Cap(_3Segment))));
        }
    }

    stSortedSet_destruct(segmentSet);

    return maximalHaplotypePaths;
}

stHash *buildSegmentToContigPathHash(stList *maximalHaplotypePaths) {
    stHash *segmentToMaximalHaplotypePathHash = stHash_construct();
    for (int32_t i = 0; i < stList_length(maximalHaplotypePaths); i++) {
        stList *maximalHaplotypePath = stList_get(maximalHaplotypePaths, i);
        assert(stList_length(maximalHaplotypePath) > 0);
        for (int32_t j = 0; j < stList_length(maximalHaplotypePath); j++) {
            Segment *segment = stList_get(maximalHaplotypePath, j);
            assert(stHash_search(segmentToMaximalHaplotypePathHash, segment)
                    == NULL);
            assert(stHash_search(segmentToMaximalHaplotypePathHash,
                    segment_getReverse(segment)) == NULL);
            stHash_insert(segmentToMaximalHaplotypePathHash, segment,
                    maximalHaplotypePath);
        }
    }
    return segmentToMaximalHaplotypePathHash;
}

int32_t contigPathLength(stList *haplotypePath) {
    int32_t k = 0;
    for (int32_t j = 0; j < stList_length(haplotypePath); j++) {
        Segment *segment = stList_get(haplotypePath, j);
        k += segment_getLength(segment);
    }
    return k;
}

stHash *buildContigPathToContigPathLengthHash(
        stList *maximalHaplotypePaths) {
    stHash *maximalHaplotypesToMaximalHaplotypePathLengths = stHash_construct();
    for (int32_t i = 0; i < stList_length(maximalHaplotypePaths); i++) {
        stList *maximalHaplotypePath = stList_get(maximalHaplotypePaths, i);
        int32_t k = contigPathLength(maximalHaplotypePath);
        stHash_insert(maximalHaplotypesToMaximalHaplotypePathLengths,
                maximalHaplotypePath, stIntTuple_construct(1, k));
    }
    return maximalHaplotypesToMaximalHaplotypePathLengths;
}

