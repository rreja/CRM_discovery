#include <stdio.h>
#include "matrix.h"
#include "alphabet.h"
#include "config.h"
#include "mcast-html.h"
#include "mhmmscan.h"
#include "motif.h"
#include "projrel.h"

const char *pick_color(int index) {
    const char *colors[] = {
      "aqua",
      "blue",
      "red",
      "fuchsia",
      "yellow",
      "lime",
      "teal",
      "#444444",
      "green",
      "silver",
      "purple",
      "olive",
      "navy",
      "maroon",
      "white"
    };

    return colors[index % 15];
}

/******************************************************************************
 * This function expects as input a string contain a sequence of DNA bases [ACGT]
 * and returns a string with HTML formatting coloring each of the bases.
 *
 * The returned string must be freed by the caller.
 * Characters not in [ACGT] will be colored black.
 *****************************************************************************/
char *color_dna_sequence(char *sequence) {

  // Each base will be converted to a longer string
  // Figure out maximum memory required
  int seq_len = strlen(sequence);
  int multiplier = strlen("<span style=\"color:orange\">c</span>");
  int mem_needed = sizeof(char) * (seq_len * multiplier + 1);
  char *colored_sequence = mymalloc(mem_needed);
  colored_sequence[0] = '\0';
  const char *color;
  while (*sequence) {
    switch (toupper(*sequence)) {
      case 'A':
        color = "red";
        break;
      case 'C':
        color = "blue";
        break;
      case 'G':
        color = "orange";
        break;
      case 'T':
        color = "green";
        break;
      default:
        color = "black";
        break;
    }
    strcat(colored_sequence, "<span style=\"color:");
    strcat(colored_sequence, color);
    strcat(colored_sequence, "\">");
    strncat(colored_sequence, sequence, 1);
    strcat(colored_sequence, "</span>");
    ++sequence;
  }
  return colored_sequence;
}

/******************************************************************************
  Print JavaScript code defining an array of motifs 
  and their best possible matches.
 *****************************************************************************/
void mcast_print_motif_array(FILE *output, MOTIF_T *motifs, int num_motifs) {
   int i;
   fputs("\n", output);
   for (i = 0; i < num_motifs; i++) {
     MOTIF_T *motif = motif_at(motifs, i);
     MOTIF_T *rc_motif = NULL;
     char *motif_id = get_motif_id(motif);
     char *rc_motif_id = NULL;
     if (i < (num_motifs - 1)) {
       rc_motif = motif_at(motifs, i + 1);
       rc_motif_id = get_motif_id(rc_motif);
     }
     char *best_possible_match = get_best_possible_match(motif);
     if (rc_motif_id && strcmp(motif_id, rc_motif_id) == 0) {
       ++i; // Pair of identiical motif ids indicate forward/reverse pair.
       char *best_possible_rc_match = get_best_possible_match(rc_motif);
       fprintf(
           output,
          "      motifs[\"%s\"] = new Motif(\"%s\", \"nucleotide\", \"%s\", \"%s\");\n",
          motif_id,
          motif_id,
          best_possible_match,
          best_possible_rc_match
       );
       myfree(best_possible_rc_match);
     }
     else {
       fprintf(
           output,
          "      motifs[\"%s\"] = new Motif(\"%s\", \"nucleotide\", \"%s\", \"%s\");\n",
          motif_id,
          motif_id,
          best_possible_match,
          "" 
       );
     }
     myfree(best_possible_match);
   }
};

void mcast_print_sequence_filename(FILE *output, MHMMSCAN_OPTIONS_T *options) {
  fputs(options->seq_filename, output);
};

void mcast_print_num_sequences(FILE *output, int num_seqs) {
  fprintf(output, "%d", num_seqs);
};

void mcast_print_num_residues(FILE *output, long num_residues) {
  fprintf(output, "%ld", num_residues);
};

void mcast_print_pattern_filename(FILE *output, MHMMSCAN_OPTIONS_T *options) {
  fputs(options->motif_filename, output);
};

void mcast_print_motif_list(FILE * output, MOTIF_T* motifs, int num_motifs) {
  fputs("\n", output);
  int i;
  for (i = 0; i < num_motifs; i++) {
    MOTIF_T *motif = motif_at(motifs, i);
    MOTIF_T *rc_motif = NULL;
    char *motif_id = get_motif_id(motif);
    int width = get_motif_length(motif);
    char *rc_motif_id = NULL;
    if (i < (num_motifs - 1)) {
      rc_motif = motif_at(motifs, i + 1);
      rc_motif_id = get_motif_id(rc_motif);
    }
    char *best_possible_match = get_best_possible_match(motif);
    char *colored_best_possible_match = color_dna_sequence(best_possible_match);
    char *best_possible_rc_match = NULL;
    char *colored_best_possible_rc_match = NULL;
    if (rc_motif_id && strcmp(motif_id, rc_motif_id) == 0) {
      ++i; // Pair of identiical motif ids indicate forward/reverse pair.
      best_possible_rc_match = get_best_possible_match(rc_motif);
      colored_best_possible_rc_match = color_dna_sequence(best_possible_rc_match);
    }
    const char *indent = "               ";
    fprintf(output, "%s<tr>\n", indent);
    fprintf(output, "%s<td>%s</td>\n", indent, motif_id);
    fprintf(output, "%s<td>%d</td>\n", indent, width);
    fprintf(output, "%s<td class=\"sequence\">%s</td>\n", indent, colored_best_possible_match);
    fprintf(output, "%s<td class=\"sequence\">%s</td>\n", indent, colored_best_possible_rc_match);
    fprintf(output, "%s</tr>\n", indent);
    myfree(best_possible_match);
    myfree(best_possible_rc_match);
    myfree(colored_best_possible_match);
    myfree(colored_best_possible_rc_match);
  }
};
 
void mcast_print_version(FILE *output) {
  fprintf(output, "%s %s", VERSION, REVISION);
};

void mcast_print_release_date(FILE *output) {
  fputs(ARCHIVE_DATE, output);
};

void mcast_print_command_line(FILE *output, MHMMSCAN_OPTIONS_T *options) {
  fputs(options->command_line, output);
};

void mcast_print_duration(FILE *output, double duration) {
  fprintf(output, "Calculation took %g seconds.", duration);
};

void mcast_print_bg_filename(FILE *output, MHMMSCAN_OPTIONS_T *options) {
    if (options->bg_filename == NULL) {
      fputs("non-redundant database", output);
    }
    else if (strcmp("--nrdb--", options->bg_filename) == 0) {
      fputs("non-redundant database", output);
    }
    else if (strcmp("--motif--", options->bg_filename) == 0) {
      fputs(options->motif_filename, output);
    }
    else if (strcmp("--uniform--", options->bg_filename) == 0) {
      fputs("uniform background", output);
    }
    else {
      fputs(options->bg_filename, output);
    }
};

void mcast_print_bg_freqs(
  FILE *output,
  ARRAY_T *bgfreqs,
  MHMMSCAN_OPTIONS_T *options
) {
  int asize = alph_size(options->alphabet, ALPH_SIZE);
  int i;
  for (i = 0; i < asize; i++) {
    if (i % 9 == 0) {
      fputc('\n', output);
    }
    fprintf(
      output,
      "%c: %1.3f ",
      alph_char(options->alphabet, i),
      get_array_item(i, bgfreqs)
    );
  }
};


