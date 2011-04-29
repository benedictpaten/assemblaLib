/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactus.h"
#include "sonLib.h"
#include "contigPaths.h"
#include "adjacencyTraversal.h"
#include "adjacencyClassification.h"

static stHash *getMaximalScaffoldPathLengthsP(stList *haplotypePaths, stHash *haplotypePathToScaffoldPathHash,
        stList *eventStrings, CapCodeParameters *capCodeParameters) {
    stHash *haplotypeToMaximalHaplotypeLengthHash = buildContigPathToContigPathLengthHash(haplotypePaths);
    stHash *segmentToMaximalHaplotypePathHash = buildSegmentToContigPathHash(haplotypePaths);
    for (int32_t i = 0; i < stList_length(haplotypePaths); i++) {
        stSortedSet *bucket = stSortedSet_construct();
        stHash_insert(haplotypePathToScaffoldPathHash, stList_get(haplotypePaths, i), bucket);
        stSortedSet_insert(bucket, stList_get(haplotypePaths, i));
    }
    for (int32_t i = 0; i < stList_length(haplotypePaths); i++) {
        stList *haplotypePath = stList_get(haplotypePaths, i);
        assert(stList_length(haplotypePath) > 0);
        Segment *_5Segment = stList_get(haplotypePath, 0);
        if (!segment_getStrand(_5Segment)) {
            _5Segment = segment_getReverse(stList_get(haplotypePath, stList_length(haplotypePath) - 1));
        }
        assert(segment_getStrand(_5Segment));
        assert(!trueAdjacency(segment_get5Cap(_5Segment), eventStrings));
        int32_t insertLength;
        int32_t deleteLength;
        enum CapCode _5CapCode = getCapCode(segment_get5Cap(_5Segment), eventStrings, &insertLength, &deleteLength, capCodeParameters);
        if (_5CapCode == SCAFFOLD_GAP || _5CapCode == AMBIGUITY_GAP) {
            assert(stHash_search(haplotypeToMaximalHaplotypeLengthHash, haplotypePath) != NULL);
            int32_t j = stIntTuple_getPosition(stHash_search(haplotypeToMaximalHaplotypeLengthHash, haplotypePath), 0);
            Segment *adjacentSegment = getAdjacentCapsSegment(segment_get5Cap(_5Segment));
            assert(adjacentSegment != NULL);
            while (!hasCapInEvents(cap_getEnd(segment_get5Cap(adjacentSegment)), eventStrings)) { //is not a haplotype end
                adjacentSegment = getAdjacentCapsSegment(segment_get5Cap(adjacentSegment));
                assert(adjacentSegment != NULL);
            }
            assert(adjacentSegment != NULL);
            assert(hasCapInEvents(cap_getEnd(segment_get5Cap(adjacentSegment)), eventStrings)); //is a haplotype end
            stList *adjacentHaplotypePath = stHash_search(segmentToMaximalHaplotypePathHash, adjacentSegment);
            if (adjacentHaplotypePath == NULL) {
                adjacentHaplotypePath = stHash_search(segmentToMaximalHaplotypePathHash, segment_getReverse(
                        adjacentSegment));
            }
            assert(adjacentHaplotypePath != NULL);
            assert(adjacentHaplotypePath != haplotypePath);
            assert(stHash_search(haplotypeToMaximalHaplotypeLengthHash, adjacentHaplotypePath) != NULL);
            int32_t k = stIntTuple_getPosition(stHash_search(haplotypeToMaximalHaplotypeLengthHash, adjacentHaplotypePath), 0);

            //Now merge the buckets and make new int tuples..
            stSortedSet *bucket1 = stHash_search(haplotypePathToScaffoldPathHash, haplotypePath);
            stSortedSet *bucket2 = stHash_search(haplotypePathToScaffoldPathHash, adjacentHaplotypePath);
            assert(bucket1 != NULL);
            assert(bucket2 != NULL);
            assert(bucket1 != bucket2);
            stSortedSet *bucket3 = stSortedSet_getUnion(bucket1, bucket2);
            stSortedSetIterator *bucketIt = stSortedSet_getIterator(bucket3);
            stList *l;
            while ((l = stSortedSet_getNext(bucketIt)) != NULL) {
                //Do the bucket first
                assert(stHash_search(haplotypePathToScaffoldPathHash, l) == bucket1 || stHash_search(haplotypePathToScaffoldPathHash, l) == bucket2);
                stHash_remove(haplotypePathToScaffoldPathHash, l);
                stHash_insert(haplotypePathToScaffoldPathHash, l, bucket3);
                //Now the length
                stIntTuple *m = stHash_remove(haplotypeToMaximalHaplotypeLengthHash, l);
                assert(m != NULL);
                assert(stIntTuple_getPosition(m, 0) == j || stIntTuple_getPosition(m, 0) == k);
                stHash_insert(haplotypeToMaximalHaplotypeLengthHash, l, stIntTuple_construct(1, j + k));
                stIntTuple_destruct(m);
            }
            assert(stHash_search(haplotypePathToScaffoldPathHash, haplotypePath) == bucket3);
            assert(stHash_search(haplotypePathToScaffoldPathHash, adjacentHaplotypePath) == bucket3);
            stSortedSet_destructIterator(bucketIt);
        }
    }
    stHash_destruct(segmentToMaximalHaplotypePathHash);
    return haplotypeToMaximalHaplotypeLengthHash;
}

static void debugMaximalHaplotypePathLengthsP(Cap *cap, stList *haplotypePath,
        stHash *haplotypePathToScaffoldPathHash, stHash *haplotypeToMaximalHaplotypeLengthHash,
        stHash *segmentToMaximalHaplotypePathHash, stList *eventStrings, CapCodeParameters *capCodeParameters, bool capDir) {
    int32_t insertLength;
    int32_t deleteLength;
    enum CapCode capCode = getCapCode(cap, eventStrings, &insertLength, &deleteLength, capCodeParameters);
    if (capCode == SCAFFOLD_GAP || capCode == AMBIGUITY_GAP) {
        Segment *adjacentSegment = getAdjacentCapsSegment(cap);
        assert(adjacentSegment != NULL);
        while (!hasCapInEvents(cap_getEnd(capDir ? segment_get5Cap(adjacentSegment) : segment_get3Cap(adjacentSegment)), eventStrings)) {
            adjacentSegment = getAdjacentCapsSegment(capDir ? segment_get5Cap(adjacentSegment) : segment_get3Cap(adjacentSegment));
            assert(adjacentSegment != NULL);
        }
        assert(adjacentSegment != NULL);
        assert(hasCapInEvents(cap_getEnd(segment_get5Cap(adjacentSegment)), eventStrings)); //isHaplotypeEnd(cap_getEnd(segment_get5Cap(adjacentSegment))));
        stIntTuple *j = stHash_search(haplotypeToMaximalHaplotypeLengthHash, haplotypePath);
        assert(j != NULL);
        stList *adjacentHaplotypePath = stHash_search(segmentToMaximalHaplotypePathHash, adjacentSegment);
        if (adjacentHaplotypePath == NULL) {
            adjacentHaplotypePath = stHash_search(segmentToMaximalHaplotypePathHash,
                    segment_getReverse(adjacentSegment));
        }
        assert(adjacentHaplotypePath != NULL);
        assert(adjacentHaplotypePath != haplotypePath);
        stIntTuple *k = stHash_search(haplotypeToMaximalHaplotypeLengthHash, adjacentHaplotypePath);
        assert(k != NULL);
        assert(stIntTuple_getPosition(j, 0) == stIntTuple_getPosition(k, 0));
        assert(stHash_search(haplotypePathToScaffoldPathHash, haplotypePath) ==
                stHash_search(haplotypePathToScaffoldPathHash, adjacentHaplotypePath));
    }
}

static void debugMaximalHaplotypePathLengths(stList *haplotypePaths, stHash *haplotypePathToScaffoldPathHash,
        stHash *haplotypeToMaximalHaplotypeLengthHash, stList *eventStrings, CapCodeParameters *capCodeParameters) {
    stHash *segmentToMaximalHaplotypePathHash = buildSegmentToContigPathHash(haplotypePaths);
    for (int32_t i = 0; i < stList_length(haplotypePaths); i++) {
        stList *haplotypePath = stList_get(haplotypePaths, i);
        assert(stList_length(haplotypePath) > 0);
        //Traversing from 5' end..
        Segment *_5Segment = stList_get(haplotypePath, 0);
        Segment *_3Segment = stList_get(haplotypePath, stList_length(haplotypePath) - 1);
        assert(segment_getStrand(_5Segment) == segment_getStrand(_3Segment));
        if (!segment_getStrand(_5Segment)) {
            Segment *j = _5Segment;
            _5Segment = segment_getReverse(_3Segment);
            _3Segment = segment_getReverse(j);
        }
        assert(segment_getStrand(_5Segment));
        assert(segment_getStrand(_3Segment));
        Cap *_5Cap = segment_get5Cap(_5Segment);
        Cap *_3Cap = segment_get3Cap(_3Segment);
        assert(!trueAdjacency(_5Cap, eventStrings));
        assert(!trueAdjacency(_3Cap, eventStrings));
        debugMaximalHaplotypePathLengthsP(_5Cap, haplotypePath,
                haplotypePathToScaffoldPathHash, haplotypeToMaximalHaplotypeLengthHash,
                segmentToMaximalHaplotypePathHash, eventStrings, capCodeParameters, 1);
        debugMaximalHaplotypePathLengthsP(_3Cap, haplotypePath,
                haplotypePathToScaffoldPathHash, haplotypeToMaximalHaplotypeLengthHash,
                segmentToMaximalHaplotypePathHash, eventStrings, capCodeParameters, 0);
    }
    stHash_destruct(segmentToMaximalHaplotypePathHash);
}

stHash *getMaximalScaffoldPathLengths(stList *haplotypePaths, stList *eventStrings, CapCodeParameters *capCodeParameters) {
    stHash *haplotypePathToScaffoldPathHash = stHash_construct();
    stHash *i = getMaximalScaffoldPathLengthsP(haplotypePaths, haplotypePathToScaffoldPathHash, eventStrings, capCodeParameters);
    debugMaximalHaplotypePathLengths(haplotypePaths, haplotypePathToScaffoldPathHash, i, eventStrings, capCodeParameters);
    stHash_destruct(haplotypePathToScaffoldPathHash);
    return i;
}

stHash *getScaffoldPaths(stList *haplotypePaths, stList *eventStrings, CapCodeParameters *capCodeParameters) {
    stHash *haplotypePathToScaffoldPathHash = stHash_construct();
    stHash *i = getMaximalScaffoldPathLengthsP(haplotypePaths, haplotypePathToScaffoldPathHash, eventStrings, capCodeParameters);
    debugMaximalHaplotypePathLengths(haplotypePaths, haplotypePathToScaffoldPathHash, i, eventStrings, capCodeParameters);
    stHash_destruct(i);
    return haplotypePathToScaffoldPathHash;
}
