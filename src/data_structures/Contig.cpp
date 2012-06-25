/*
 * Contig.cpp
 *
 *  Created on: Jan 25, 2012
 *      Author: vezzi
 */

#include "Contig.h"


Contig::Contig() {

}

Contig::Contig(unsigned int contigLength, unsigned int libraries) {
	this->contigLength = contigLength;
	this->libraries = libraries;
	this->readCoverage = new unsigned int*[libraries+1];
	this->spanCoverage = new unsigned int*[libraries+1];
	this->maxInserts   = new unsigned int[libraries];
	this->minInserts   = new unsigned int[libraries];
	for(unsigned int i=0; i < libraries; i ++) {
		this->readCoverage[i] = new unsigned int[contigLength];
		this->spanCoverage[i] = new unsigned int[contigLength];
		this->minInserts[i] = 0;
		this->maxInserts[i] = 10000000;
	}
	this->readCoverage[libraries] = new unsigned int[contigLength];
	this->spanCoverage[libraries] = new unsigned int[contigLength];

	// reset memory location WRITE 0 everywhere
	for(unsigned int i = 0; i < libraries + 1 ; i ++) {
		for(unsigned int j = 0; j < contigLength ; j++) {
			this->readCoverage[i][j] = 0;
			this->spanCoverage[i][j] = 0;
		}
	}



}

Contig::~Contig() {
	if(readCoverage != NULL) {
		for(unsigned int i=0 ; i < libraries ; i++) {
			delete [] this->readCoverage[i];
		}
		delete [] this->readCoverage;
	}
}


void Contig::setLibraryLimits(unsigned int minInsert, unsigned int maxInsert, unsigned int library) {
	if(library >= 0 and library < this->libraries) {
		this->minInserts[library] = minInsert;
		this->maxInserts[library] = maxInsert;
	} else {
		cout << "in Contig::setLibrariesLimits passed a library identifier larger than the number of libraries\n";
		exit(1);
	}
}


void Contig::updateReadCov(unsigned int start, unsigned int end, unsigned int library) {
	if(start < 0) {
		cout << "hoops, start less than 0 when updating CONTIG\n";
		start = 0;
	}
	if(end > this->contigLength) {
		end = this->contigLength;
	}
	// now update
	for(unsigned int i = start; i< end ; i ++ ) {
		this->readCoverage[library][i]++;
		this->readCoverage[this->libraries][i]++;
	}

}

void Contig::updateSpanCov(unsigned int start, unsigned int end, unsigned int library) {
	if(start < 0) {
		cout << "hoops, start less than 0 when updating CONTIG\n";
		start = 0;
	}
	if(end > this->contigLength) {
		end = this->contigLength;
	}
	// now update
	for(unsigned int i = start; i< end ; i ++ ) {
		this->spanCoverage[library][i]++;
		this->spanCoverage[this->libraries][i]++;
	}

}


void Contig::updateContig(bam1_t* b, unsigned int library) {
	const bam1_core_t* core =  &b->core;
	uint32_t* cigar = bam1_cigar(b);
	int32_t alignmentLength = 0;
	int32_t startRead=0;
	int32_t endRead=0;
	int32_t startPaired=0;
	int32_t startInsert=0;
	int32_t endInsert=0;
	uint32_t iSize=0;

	uint32_t minInsert = this->minInserts[library];
	uint32_t maxInsert = this->maxInserts[library];

	if(!(core->flag&BAM_FUNMAP) && !(core->flag&BAM_FDUP) && !(core->flag&BAM_FSECONDARY) && !(core->flag&BAM_FQCFAIL)) { // if read has been mapped and it is not a DUPLICATE or a SECONDARY alignment
		alignmentLength = bam_cigar2qlen(core,cigar);
		startRead = core->pos; // start position on the contig
		endRead = startRead + alignmentLength ; // position where reads ends
		updateReadCov(startRead, endRead, library); //  update read coverage

		iSize = abs(core->isize);

		if ((core->flag&BAM_FREAD1) //First in pair
				&& !(core->flag&BAM_FMUNMAP) /*Mate is also mapped!*/
				&& (core->tid == core->mtid) /*Mate on the same chromosome*/
		) {
			startPaired = core->mpos;
			if(startRead < startPaired) {
				iSize = (startPaired + core->l_qseq -1) - startRead; // insert size, I consider both reads of the same length
				startInsert = startRead;
				endInsert = startRead + iSize;
				if (minInsert <= iSize && iSize <= maxInsert) { // DO NOT CHECK ORIENTATION
					updateSpanCov(startInsert,endInsert, library); // update spanning coverage
				}

				if(!(core->flag&BAM_FREVERSE) && (core->flag&BAM_FMREVERSE) ) { //
					//here reads are correctly oriented
					if (minInsert <= iSize && iSize <= maxInsert) { //this is a right insert
						//updateSpanCov(startRead, endRead, library); // update good read coverage
					}
				}
			} else {
				iSize = (startRead + alignmentLength - 1) - startPaired;
				startInsert = startPaired;
				endInsert = startInsert + iSize;

				if (minInsert <= iSize && iSize <= maxInsert) { //
					updateSpanCov(startInsert,endInsert, library); // update spanning coverage
				}

				if((core->flag&BAM_FREVERSE) && !(core->flag&BAM_FMREVERSE) ) { //
					//here reads are correctly oriented
					if (minInsert <= iSize && iSize <= maxInsert) { //this is a right insert
						//updateSpanCov(startRead, endRead, cmCov); // update good read coverage
					}
				}
			}
		} // other case are not considered here

	}
}




