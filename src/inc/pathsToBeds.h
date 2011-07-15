/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef PATHS_TO_BEDS_H_
#define PATHS_TO_BEDS_H_

#include "cactus.h"
#include "sonLib.h"
#include "adjacencyClassification.h"
#include "contigPaths.h"
#include "scaffoldPaths.h"

typedef struct _sequenceInterval {
        int32_t start;
        int32_t end;
        char *sequenceName;
} SequenceInterval;

SequenceInterval *sequenceInterval_construct(int32_t start, int32_t end,
        const char *sequenceName);

void sequenceInterval_destruct(SequenceInterval *sequenceInterval);

stList *getContigPathIntervals(Flower *flower, stList *contigPaths, const char *chosenEventString, stList *referenceEventStrings);

stList *getScaffoldPathIntervals(Flower *flower, const char *chosenEventString, stList *referenceEventStrings, stList *contaminationEventStrings, CapCodeParameters *capCodeParameters);

stList *getSplitContigPathIntervals(Flower *flower, stList *contigPaths, const char *chosenEventString,
        stList *eventStrings);

#endif
