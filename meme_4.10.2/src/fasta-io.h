/****************************************************************************
 * FILE: fasta-io.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 9/14/98
 * PROJECT: MHMM
 * DESCRIPTION: Read sequences from a FASTA file into memory.
 * COPYRIGHT: 1998-2008 WSN
 ****************************************************************************/
#ifndef FASTA_IO_H
#define FASTA_IO_H
#include "seq.h"
#include "seq-reader-from-fasta.h"

#define MAX_SEQ 250000000   // Default maximum sequence length

/****************************************************************************
 * Read one sequence from a file in Fasta format.
 *
 * Return: Was a sequence successfully read?
 ****************************************************************************/
BOOLEAN_T read_one_fasta
  (ALPH_T    alph,
   FILE*     fasta_file,
   int       max_seq,
   SEQ_T**   sequence);

/****************************************************************************
 * Read up to N letters of one sequence from a file in Fasta format.
 * On subsequent calls, read either
 *  - the next N letters from the same sequence, or 
 *  - the first N letters from the next sequence (if the whole
 *    sequence was read last time).
 * In the second case, the new segment is appended to the end of the
 * existing sequence segment.
 *
 * Return: Was a sequence segment successfully read?
 ****************************************************************************/
BOOLEAN_T read_one_fasta_segment(
   ALPH_T alph,
   FILE* fasta_file,
   int max_chars,
   SEQ_T** sequence
);

/****************************************************************************
 * Read up to max_chars letters of one sequence from a DATA_BLOCK_T readder
 * and copy them in to the raw sequence in the SEQ_T object starting at the
 * given buffer offset. 
 ****************************************************************************/
void read_one_fasta_segment_from_reader(
  DATA_BLOCK_READER_T *fasta_reader,
  size_t max_size,
  size_t buffer_offset,
  SEQ_T* sequence
);


/****************************************************************************
 * Read all the sequences from a file in Fasta format.
 ****************************************************************************/
void read_many_fastas
  (ALPH_T     alph,
   FILE*      fasta_file,
   int        max_seq, // Maximum sequence length.
   int*       num_seqs,
   SEQ_T***   sequences);

/****************************************************************************
 *  Create a zero order background from the passed sequences
 ****************************************************************************/
ARRAY_T* calc_bg_from_fastas
  (ALPH_T alph, 
   int num_seqs, 
   SEQ_T** sequences);

/****************************************************************************
 * Print a single sequence in FASTA format.
 ****************************************************************************/
void print_fasta
  (SEQ_T*    sequence,
   FILE*     outfile);

/****************************************************************************
 *  Create a sequence object from the first sequence in a FASTA file.
 ****************************************************************************/
SEQ_T* read_sequence_from_file
(ALPH_T alph,
 char* filename);

/****************************************************************************
 *  Create a zero order background from the sequences in a FASTA file.
 ****************************************************************************/
ARRAY_T* calc_bg_from_file
  (ALPH_T alph, 
   char* filename, 
   BOOLEAN_T allow_stdin);

#endif
