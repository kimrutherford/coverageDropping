/*
 * Contig.h
 *
 *  Created on: Jan 25, 2012
 *      Author: vezzi
 */

#ifndef CONTIG_H_
#define CONTIG_H_

using namespace std;

#include "samtools/sam.h"
#include "options/Options.h"
#include "common.h"



class Contig{
	unsigned int contigLength;
	unsigned int libraries;
	unsigned int *minInserts;
	unsigned int *maxInserts;
	unsigned int *meanInserts;
	unsigned int *edgeCutoff;
	unsigned int *dropCutoff;
	unsigned int **spanCoverage;
	unsigned int **readCoverage;

	unsigned long int *totalReadLength;
	unsigned long int *totalSpanLength;

	float *readCov;
	float *spanCov;

	void updateReadCov(unsigned int start, unsigned int end, unsigned int library);
	void updateSpanCov(unsigned int start, unsigned int end, unsigned int library);


public:
	Contig();
	Contig(unsigned int contigLength, unsigned int libraries);
	~Contig();

	void setUpTotal();
	void setLibraryLimits(unsigned int library, unsigned int minInsert, unsigned int maxInsert, unsigned int meanInsert,
			 unsigned int edgeCutoff,  unsigned int dropCutoff); // set min and max of liobrary library

	void updateContig(const bam1_t* b, unsigned int library); // given an alignment and its library it updates the contig situation

	void computeContigStats();

	void computeCoverageDrops();

	void computeCoverageDrops(string ctgID, unsigned int lib, ofstream &output);

	void print();
	void printStats();
	void printLibrariesThresholds();

};





#endif /* CONTIG_H_ */
