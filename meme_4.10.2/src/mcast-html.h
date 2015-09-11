#ifndef MCAST_HTML_H
#define MCAST_HTML_H
#include <stdio.h>
#include "mhmmscan.h"
#include "array.h"

const char *pick_color(int index);

void mcast_print_sequence_filename(FILE *output, MHMMSCAN_OPTIONS_T *options);
void mcast_print_num_sequences(FILE *output, int num_seqs);
void mcast_print_num_residues(FILE *output, long num_residues);
void mcast_print_pattern_filename(FILE *output, MHMMSCAN_OPTIONS_T *options);
void mcast_print_motif_array(FILE *output, MOTIF_T *motifs, int num_motifs);
void mcast_print_motif_list(FILE * output, MOTIF_T* motifs, int num_motifs);
void mcast_print_version(FILE *output);
void mcast_print_release_date(FILE *output);
void mcast_print_command_line(FILE *output, MHMMSCAN_OPTIONS_T *options);
void mcast_print_duration(FILE *output, double duration);
void mcast_print_bg_filename(FILE *output, MHMMSCAN_OPTIONS_T *options);
void mcast_print_bg_freqs( FILE *output, ARRAY_T *bgfreqs, MHMMSCAN_OPTIONS_T *options);
#endif
