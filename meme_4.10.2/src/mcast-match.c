#include <assert.h>
#include <string.h>
#include "mcast-html-string.h"
#include "mcast-html.h"
#include "mcast-match.h"
#include "object-list.h"
#include "utils.h"
#include "cisml.h"

struct motif_hit_t {
  char *motif_id;
  int motif_index;
  char *seq;
  size_t start;
  size_t stop;
  char strand;
  double pvalue;
};

struct mcast_match_t {
  int cluster_id;
  OBJECT_LIST_T *motif_hits;
  char *seq_name;
  char *sequence;
  int gc_bin;
  double score;
  double gc;
  double evalue;
  double pvalue;
  double qvalue;
  size_t start;
  size_t stop;
};

/***********************************************************************
 * allocate_mcast_match
 *
 * This function generates a new mcast match object.
 ***********************************************************************/
MCAST_MATCH_T *allocate_mcast_match() {
  static int cluster_id = 0;
  MCAST_MATCH_T *mcast_match = mm_malloc(sizeof(MCAST_MATCH_T));
  ++cluster_id;
  mcast_match->cluster_id = cluster_id;

  mcast_match->motif_hits
    = new_object_list(NULL, NULL, NULL, free_motif_hit);
  mcast_match->seq_name = NULL;
  mcast_match->sequence = NULL;
  mcast_match->gc_bin = -1;
  mcast_match->score = -1.0;
  mcast_match->gc = NaN();
  mcast_match->evalue = NaN();
  mcast_match->pvalue = NaN();
  mcast_match->qvalue = NaN();
  mcast_match->start = -1;
  mcast_match->stop = -1;

  return mcast_match;
}

/***********************************************************************
 * copy_mcast_match
 *
 * This function deep copies an mcast match object.
 *
 ***********************************************************************/
void *copy_mcast_match(void *m) {
  MCAST_MATCH_T *mcast_match = m;
  MCAST_MATCH_T *new_mcast_match = mm_malloc(sizeof(MCAST_MATCH_T));

  new_mcast_match->motif_hits
    = new_object_list(NULL, NULL, NULL, free_motif_hit);
  MOTIF_HIT_T *motif_hit = retrieve_next_object(mcast_match->motif_hits);
  while (motif_hit != NULL) {
    // FIXME actually do the copying
    motif_hit = retrieve_next_object(mcast_match->motif_hits);
  }

  new_mcast_match->seq_name = strdup(mcast_match->seq_name);
  new_mcast_match->sequence = strdup(mcast_match->sequence);
  new_mcast_match->gc_bin = mcast_match->gc_bin;
  new_mcast_match->score = mcast_match->score;
  new_mcast_match->gc = mcast_match->gc;
  new_mcast_match->evalue = mcast_match->evalue;
  new_mcast_match->pvalue = mcast_match->pvalue;
  new_mcast_match->qvalue = mcast_match->qvalue;
  new_mcast_match->start = mcast_match->start;
  new_mcast_match->stop = mcast_match->stop;

  return new_mcast_match;
}

/***********************************************************************
 * free_mcast_match
 *
 * This function frees the memory associated with a mcast_matc object.
 ***********************************************************************/
void free_mcast_match(void* m) {
  MCAST_MATCH_T *mcast_match = m;
  free_object_list(mcast_match->motif_hits);
  myfree(mcast_match->seq_name);
  myfree(mcast_match->sequence);
  myfree(mcast_match);
}

/***********************************************************************
 * compare_mcast_matches
 *
 * This function compares two mcast match objects.
 * mcast match objects are ordered by pvalue from lower to higher.
 * This function returns 
 *    1 if m1->pvalue > m2->pvalue,
 *    0 if m1->pvalue = m2->pvalue,
 *    -1 if m1->pvalue < m2->pvalue,
 *
 ***********************************************************************/
int compare_mcast_matches(void *m1, void *m2) {
  MCAST_MATCH_T *match1 = m1;
  MCAST_MATCH_T *match2 = m2;

  if (match1->pvalue > match2->pvalue) {
    return -1;
  }
  else if (match1->pvalue == match2->pvalue) {
    return 0;
  }
  else {
    return 1;
  }
}

/***********************************************************************
 * compare_mcast_match_pvalues
 *
 * This function compares the p-values of two MCAST matches
 * and returns -1, 0, or 1 as the p-value of the 1st argument is
 * less than, equal to, or greater than the p-value of the 2nd argument.
 *
 ***********************************************************************/
int compare_mcast_match_pvalues(const void *param1, const void *param2) {
  MCAST_MATCH_T *match1 = *((MCAST_MATCH_T **) param1);
  MCAST_MATCH_T *match2 = *((MCAST_MATCH_T **) param2);

  if (match1->pvalue > match2->pvalue) {
    return 1;
  }
  else if (match1->pvalue < match2->pvalue) {
    return -1;
  }
  else {
    return 0;
  }
}

/***********************************************************************
 * rev_compare_mcast_match_pvalues
 *
 * This function compares the p-values of two MCAST matches
 * and returns 1, 0, or -1 as the p-value of the 1st argument is
 * less than, equal to, or greater than the p-value of the 2nd argument.
 *
 ***********************************************************************/
int rev_compare_mcast_match_pvalues(const void *param1, const void *param2) {
  MCAST_MATCH_T *match1 = *((MCAST_MATCH_T **) param1);
  MCAST_MATCH_T *match2 = *((MCAST_MATCH_T **) param2);

  if (match1->pvalue > match2->pvalue) {
    return -1;
  }
  else if (match1->pvalue < match2->pvalue) {
    return 1;
  }
  else {
    return 0;
  }
}

/***********************************************************************
 * get_mcast_match_cluster_id
 *
 * This function returns the cluster id from an MCAST match object.
 *
 ***********************************************************************/
int get_mcast_match_cluster_id(MCAST_MATCH_T *match) {
  return match->cluster_id;
}

/***********************************************************************
 * set_mcast_match_seq_name
 *
 * This function sets the seq_name for an MCAST match object.
 *
 * The string will be duplicated and freed when the mcatch object 
 * is destroyed.
 *
 ***********************************************************************/
void set_mcast_match_seq_name(MCAST_MATCH_T *match, char *name) {
  match->seq_name = strdup(name);
}

/***********************************************************************
 * get_mcast_match_seq_name
 *
 * This function returns the seq_name from an MCAST match object.
 *
 * The caller should NOT free the returned string. It will be freed
 * when the match object is destroyed.
 *
 ***********************************************************************/
char *get_mcast_match_seq_name(MCAST_MATCH_T *match) {
  return match->seq_name;
}

/***********************************************************************
 * set_mcast_match_sequence
 *
 * This function sets the sequence for an MCAST match object.
 *
 * The string will be duplicated and freed when the mcatch object 
 * is destroyed.
 *
 ***********************************************************************/
void set_mcast_match_sequence(MCAST_MATCH_T *match, char *sequence) {
  match->sequence = strdup(sequence);
}

/***********************************************************************
 * get_mcast_match_sequence
 *
 * This function returns the sequence from an MCAST match object.
 *
 * The caller should NOT free the returned string. It will be freed
 * when the match object is destroyed.
 *
 ***********************************************************************/
char *get_mcast_match_sequence(MCAST_MATCH_T *match) {
  return match->sequence;
}

/***********************************************************************
 * set_mcast_match_gc_bin
 *
 * This function sets the gc_bin for an MCAST match object.
 *
 ***********************************************************************/
void set_mcast_match_gc_bin(MCAST_MATCH_T *match, int gc_bin) {
  match->gc_bin = gc_bin;
}

/***********************************************************************
 * get_mcast_match_gc_bin
 *
 * This function gets the gc_bin for an MCAST match object.
 *
 ***********************************************************************/
int get_mcast_match_gc_bin(MCAST_MATCH_T *match) {
  return match->gc_bin;
}

/***********************************************************************
 * set_mcast_match_score
 *
 * This function sets the score for an MCAST match object.
 *
 ***********************************************************************/
void set_mcast_match_score(MCAST_MATCH_T *match, double score) {
  match->score = score;
}

/***********************************************************************
 * get_mcast_match_score
 *
 * This function returns the score from an MCAST match object.
 *
 ***********************************************************************/
double get_mcast_match_score(MCAST_MATCH_T *match) {
  return match->score;
}

/***********************************************************************
 * set_mcast_match_gc
 *
 * This function sets the gc value for an MCAST match object.
 *
 ***********************************************************************/
void set_mcast_match_gc(MCAST_MATCH_T *match, double gc) {
  match->gc = gc;
}

/***********************************************************************
 * get_mcast_match_gc
 *
 * This function returns the gc value from an MCAST match object.
 *
 ***********************************************************************/
double get_mcast_match_gc(MCAST_MATCH_T *match) {
  return match->gc;
}

/***********************************************************************
 * set_mcast_match_evalue
 *
 * This function sets the evalue for an MCAST match object.
 *
 ***********************************************************************/
void set_mcast_match_evalue(MCAST_MATCH_T *match, double evalue) {
  match->evalue = evalue;
}

/***********************************************************************
 * get_mcast_match_evalue
 *
 * This function returns the evalue from an MCAST match object.
 *
 ***********************************************************************/
double get_mcast_match_evalue(MCAST_MATCH_T *match) {
  return match->evalue;
}

/***********************************************************************
 * set_mcast_match_pvalue
 *
 * This function sets the pvalue for an MCAST match object.
 *
 ***********************************************************************/
void set_mcast_match_pvalue(MCAST_MATCH_T *match, double pvalue) {
  match->pvalue = pvalue;
}

/***********************************************************************
 * get_mcast_match_pvalue
 *
 * This function returns the pvalue from an MCAST match object.
 *
 ***********************************************************************/
double get_mcast_match_pvalue(MCAST_MATCH_T *match) {
  return match->pvalue;
}

/***********************************************************************
 * set_mcast_match_qvalue
 *
 * This function sets the qvalue for an MCAST match object.
 *
 ***********************************************************************/
void set_mcast_match_qvalue(MCAST_MATCH_T *match, double qvalue) {
  match->qvalue = qvalue;
}

/***********************************************************************
 * get_mcast_match_qvalue
 *
 * This function returns the qvalue from an MCAST match object.
 *
 ***********************************************************************/
double get_mcast_match_qvalue(MCAST_MATCH_T *match) {
  return match->qvalue;
}

/***********************************************************************
 * set_mcast_match_start
 *
 * This function sets the match start position for an MCAST match object.
 *
 ***********************************************************************/
void set_mcast_match_start(MCAST_MATCH_T *match, size_t start) {
  match->start = start;
}

/***********************************************************************
 * get_mcast_match_start
 *
 * This function returns the match start position from an MCAST match object.
 *
 ***********************************************************************/
size_t get_mcast_match_start(MCAST_MATCH_T *match) {
  return match->start;
}

/***********************************************************************
 * set_mcast_match_stop
 *
 * This function sets the match stop position for an MCAST match object.
 *
 ***********************************************************************/
void set_mcast_match_stop(MCAST_MATCH_T *match, size_t stop) {
  match->stop = stop;
}

/***********************************************************************
 * get_mcast_match_stop
 *
 * This function returns the match start position from an MCAST match object.
 *
 ***********************************************************************/
size_t get_mcast_match_stop(MCAST_MATCH_T *match) {
  return match->stop;
}

/***********************************************************************
 * allocate_motif_hit
 *
 * This function generates a new motif hit object.
 ***********************************************************************/
MOTIF_HIT_T *allocate_motif_hit(
  char *motif_id, 
  int motif_index,
  char *seq,
  char strand, 
  size_t start, 
  size_t stop, 
  double pvalue
) {
  MOTIF_HIT_T *motif_hit = mm_malloc(sizeof(MOTIF_HIT_T));

  motif_hit->motif_id = strdup(motif_id);
  motif_hit->motif_index = motif_index;
  motif_hit->seq = strdup(seq);
  motif_hit->strand = strand;
  motif_hit->start = start;
  motif_hit->stop = stop;
  motif_hit->pvalue = pvalue;

  return motif_hit;
}

/***********************************************************************
 * free_motif_hit
 *
 * This function frees the memory associated with a motif hit object.
 ***********************************************************************/
void free_motif_hit(MOTIF_HIT_T* motif_hit) {
  myfree(motif_hit->motif_id);
  myfree(motif_hit->seq);
  myfree(motif_hit);
}

/***********************************************************************
 * get_motif_hit_motif_id
 *
 * This function returns a pointer to a string containing the motif id
 * associated with a motif hit object. 
 *
 * The caller is NOT responsible for freeing the returnted string.
 * It will be freed when the motif hit is freed.
 *
 ***********************************************************************/
const char* get_motif_hit_motif_id(MOTIF_HIT_T* motif_hit) {
  return motif_hit->motif_id;
}

/***********************************************************************
 * get_motif_hit_seq
 *
 * This function returns a pointer to a string containing the motif id
 * associated with a motif hit object. 
 *
 * The caller is NOT responsible for freeing the returnted string.
 * It will be freed when the motif hit is freed.
 *
 ***********************************************************************/
const char* get_motif_hit_seq(MOTIF_HIT_T* motif_hit) {
  return motif_hit->seq;
}

/***********************************************************************
 * get_motif_hit_strand
 *
 * This function returns the strand assocaited with a motif hit object.
 *
 ***********************************************************************/
char get_motif_hit_strand(MOTIF_HIT_T* motif_hit) {
  return motif_hit->strand;
}

/***********************************************************************
 * get_motif_hit_start
 *
 * This function returns the starting position assocaited with 
 * a motif hit object.
 *
 ***********************************************************************/
size_t get_motif_hit_start(MOTIF_HIT_T* motif_hit) {
  return motif_hit->start;
}

/***********************************************************************
 * get_motif_hit_stop
 *
 * This function returns the stop position assocaited with 
 * a motif hit object.
 *
 ***********************************************************************/
size_t get_motif_hit_stop(MOTIF_HIT_T* motif_hit) {
  return motif_hit->stop;
}

/***********************************************************************
 * get_motif_hit_pvalue
 *
 * This function returns the pvalue assocaited with a motif hit object.
 *
 ***********************************************************************/
double get_motif_hit_pvalue(MOTIF_HIT_T* motif_hit) {
  return motif_hit->pvalue;
}

/***********************************************************************
 * add_mcast_match_motif_hit
 *
 * This function adds a motif_hit to an mcast_match.
 * Freeing the mcast_match will free the motif_hits that have been added.
 ***********************************************************************/
void add_mcast_match_motif_hit(
    MCAST_MATCH_T* mcast_match,
    MOTIF_HIT_T *motif_hit
) {
  store_object(motif_hit, NULL, 0.0, mcast_match->motif_hits);
}

/***********************************************************************
 * print_mcast_match
 *
 * This function prints the data from mcast_match object
 *
 ***********************************************************************/
void print_mcast_match(FILE *output, void *m) {
  MCAST_MATCH_T *mcast_match = m;
  fprintf(
    output,
    "%s\t%5.1g\t%3.1g\t%3.1g\t%zd\t%zd\n",
    mcast_match->seq_name,
    mcast_match->score,
    mcast_match->pvalue,
    mcast_match->evalue,
    mcast_match->start,
    mcast_match->stop
  );
  MOTIF_HIT_T *motif_hit = retrieve_next_object(mcast_match->motif_hits);
  while(motif_hit != NULL) {
    fprintf(
      output,
      "\t%s\t%3.1g\t%c\t%zd\t%zd\n",
      motif_hit->motif_id,
      motif_hit->pvalue,
      motif_hit->strand,
      motif_hit->start,
      motif_hit->stop
    );
    motif_hit = retrieve_next_object(mcast_match->motif_hits);
  }
}
/**********************************************************************
  print_mcast_match_as_cisml

  Print a heap of mcast_matches asCisML XML
**********************************************************************/
void print_mcast_match_as_cisml(
  FILE* out, 
  MCAST_MATCH_T *mcast_match
) {

  assert(out != NULL);
  assert(mcast_match != NULL);

  int cluster_id = get_mcast_match_cluster_id(mcast_match);
  char *seq_name = get_mcast_match_seq_name(mcast_match);
  char *match_seq = get_mcast_match_sequence(mcast_match);
  double match_score = get_mcast_match_score(mcast_match);
  double match_evalue = get_mcast_match_evalue(mcast_match);
  double match_pvalue = get_mcast_match_pvalue(mcast_match);
  double match_qvalue = get_mcast_match_qvalue(mcast_match);
  size_t match_start = get_mcast_match_start(mcast_match);
  size_t match_stop = get_mcast_match_stop(mcast_match);

  fprintf(out, "<multi-pattern-scan");
  fprintf(out, " score=\"%g\"", match_score);
  fprintf(out, " pvalue=\"%.5g\"", match_pvalue);
  fprintf(out, ">\n");
  MOTIF_HIT_T *motif_hit = retrieve_next_object(mcast_match->motif_hits);
  while(motif_hit != NULL) {
    const char *motif_id = get_motif_hit_motif_id(motif_hit);
    const char *hit_seq = get_motif_hit_seq(motif_hit);
    char strand = get_motif_hit_strand(motif_hit);
    double hit_pvalue = get_motif_hit_pvalue(motif_hit);
    size_t hit_start = get_motif_hit_start(motif_hit);
    size_t hit_stop = get_motif_hit_stop(motif_hit);
    if (strand == '-') {
      // Reverse strand, swap start and stop
      size_t tmp = hit_stop;
      hit_stop = hit_start;
      hit_start = tmp;
    }
    fprintf(
      out,
      "<pattern accession=\"%s\" name=\"%s\">\n",
      motif_id,
      motif_id
    );
    fprintf(
      out,
      "<scanned-sequence accession=\"%s\" name=\"%s\">\n",
      seq_name,
      seq_name
    );
    fprintf(
      out,
      "<matched-element start=\"%zd\" stop=\"%zd\" pvalue=\"%.5g\">\n",
      hit_start,
      hit_stop,
      hit_pvalue
    );
    fprintf(out, "<sequence>%s</sequence>\n", hit_seq);
    fputs("</matched-element>\n", out);
    fputs("</scanned-sequence>\n", out);
    fputs("</pattern>\n", out);
    motif_hit = retrieve_next_object(mcast_match->motif_hits);
  }
  fprintf(
    out, 
    "<mem:match cluster-id=\"cluster-%d\" "
    "seq-name=\"%s\" start=\"%zd\" stop=\"%zd\" "
    "evalue=\"%.5g\" qvalue=\"%.5g\""
    ">",
    cluster_id,
    seq_name,
    match_start,
    match_stop,
    match_evalue,
    match_qvalue
  );
  fprintf(out, "%s\n", match_seq);
  fputs("</mem:match>\n", out);
  fputs("</multi-pattern-scan>\n", out);
}

/**********************************************************************
  mcast_print_results_as_cisml

  Print a heap of mcast_matches as CisML XML
**********************************************************************/
void mcast_print_results_as_cisml(
  BOOLEAN_T stats_available,
  int num_matches,
  MCAST_MATCH_T **mcast_matches,
  MHMMSCAN_OPTIONS_T *options
) {

  assert(options != NULL);
  assert(mcast_matches != NULL);

  // Open the CisML XML file for output
  FILE *out = fopen(options->cisml_path, "w");
  if (!out) {
    die("Couldn't open file %s for output.\n", options->cisml_path);
  }

  print_cisml_start(out, options->program, TRUE, NULL, TRUE);

  fprintf(out, "<parameters>\n");
  fprintf(
    out,
    "<pattern-file>%s</pattern-file>\n",
    options->motif_filename 
  );
  fprintf(
    out,
    "<sequence-file>%s</sequence-file>\n",
    options->seq_filename
  );
  if (options->bg_filename != NULL) {
    fprintf(
      out,
      "<background-seq-file>%s</background-seq-file>\n",
      options->bg_filename
    );
  }
  fprintf(
    out,
    "<pattern-pvalue-cutoff>%g</pattern-pvalue-cutoff>\n",
    options->motif_pthresh
  );
  fprintf(
    out,
    "<sequence-pvalue-cutoff>%g</sequence-pvalue-cutoff>\n",
    options->p_thresh
  );
  fprintf(out, "</parameters>\n" );

  int i = 0;
  for (i = 0; i < num_matches; ++i) {
    MCAST_MATCH_T *match = mcast_matches[i];
    double evalue = get_mcast_match_evalue(match);
    double pvalue = get_mcast_match_pvalue(match);
    double qvalue = get_mcast_match_qvalue(match);
    if (!stats_available
       || (
         evalue <= options->e_thresh 
         && pvalue <= options->p_thresh 
         && qvalue <= options->q_thresh)
     ) {
      print_mcast_match_as_cisml(out, match);
    }
  }

  print_cisml_end(out);

  fclose(out);
}

/**********************************************************************
  mcast_print_results_as_gff

  Print an array of mcast_matches as gff.
**********************************************************************/
void mcast_print_results_as_gff(
  BOOLEAN_T stats_available,
  int num_matches,
  MCAST_MATCH_T **mcast_matches,
  MHMMSCAN_OPTIONS_T *options
) {
  assert(options != NULL);
  assert(mcast_matches != NULL);

  FILE *out = fopen(options->gff_path, "w");
  if (!out) {
    die("Couldn't open file %s for output.\n", options->gff_path);
  }

  // Print header line
  char *header = "##gff-version 3\n";
  fprintf(out, "%s", header);

  int i = 0;
  for (i = 0; i < num_matches; ++i) {
    MCAST_MATCH_T *match = mcast_matches[i];
    double evalue = get_mcast_match_evalue(match);
    double pvalue = get_mcast_match_pvalue(match);
    double qvalue = get_mcast_match_qvalue(match);
    if (!stats_available
       || (
         evalue <= options->e_thresh 
         && pvalue <= options->p_thresh 
         && qvalue <= options->q_thresh)
     ) {
      int id = get_mcast_match_cluster_id(match);
      char *seq_name = get_mcast_match_seq_name(match);
      char *seq = get_mcast_match_sequence(match);
      double score = MIN(1000.0, -4.342 * my_log(pvalue)); // -10 * log10(pvalue)
      size_t start = get_mcast_match_start(match);
      size_t stop = get_mcast_match_stop(match);
      fprintf(
        out,
        "%s\tMCAST\ttranscriptional_cis_regulatory_region\t%zd\t%zd\t%3.3g\t.\t.\t"
        "ID=cluster-%d;pvalue=%.3g;evalue=%3.3g;qvalue=%.3g\n",
        seq_name,
        start,
        stop,
        score,
        id,
        pvalue,
        evalue,
        qvalue
      );
    }
  }

  fclose(out);
}

void mcast_print_longest_seq_name(
  FILE *output,
  MCAST_MATCH_T **mcast_matches,
  int num_matches
) {
  int max_length = 0;
  int i = 0;
  for (i = 0; i < num_matches; ++i) {
    int length = strlen(mcast_matches[i]->seq_name);
    if (length >= max_length) {
      max_length = length;
    }
  }
  fprintf(output, "%gem", 0.8 * max_length);
}

size_t get_mcast_matches_max_match_length(
  MCAST_MATCH_T **mcast_matches,
  int num_matches
) {
  size_t max_length = 0;
  int i = 0;
  for (i = 0; i < num_matches; ++i) {
    size_t length = mcast_matches[i]->stop - mcast_matches[i]->start + 1;
    if (length >= max_length) {
      max_length = length;
    }
  }

  return max_length;
}

void mcast_print_max_match_length(
  FILE *output,
  MCAST_MATCH_T **mcast_matches,
  int num_matches
) {
  size_t max_length = get_mcast_matches_max_match_length(mcast_matches, num_matches);
  fprintf(output, "%ld", max_length);
}

// Requires that statistics are available
void count_passing_matches(
  MCAST_MATCH_T **mcast_matches,
  int num_matches,
  MHMMSCAN_OPTIONS_T *options,
  int *num_passing_matches,
  double *output_thresh
) {

  switch(options->output_thresh_type) {
    case PVALUE:
      *output_thresh = options->p_thresh;
      break;
    case EVALUE:
      *output_thresh = options->e_thresh;
      break;
    case QVALUE:
      *output_thresh = options->q_thresh;
      break;
    default:
      die("Unrecognized output threshold type.");
      break;
  }
  double value = NaN();
  int i;
  for (i = 0; i < num_matches; ++i) {
    switch(options->output_thresh_type) {
      case PVALUE:
        value = get_mcast_match_pvalue(mcast_matches[i]);
        break;
      case EVALUE:
        value = get_mcast_match_evalue(mcast_matches[i]);
        break;
      case QVALUE:
        value = get_mcast_match_qvalue(mcast_matches[i]);
        break;
      default:
        die("Unrecognized output threshold type.");
        break;
    }
    if (value > *output_thresh) { break; }
  }
  *num_passing_matches = i;


}

void mcast_print_matches(
  FILE *output,
  MCAST_MATCH_T **mcast_matches,
  int num_matches,
  BOOLEAN stats_available,
  MOTIF_T *motifs,
  int num_motifs,
  MHMMSCAN_OPTIONS_T *options
) {

  const char *output_thresh_type[] = {
    "a <i>p</i>-value",
    "an <i>E</i>-value",
    "a <i>q</i>-value"
  };

  int num_passing_matches = num_matches;
  if (stats_available) {
    double output_thresh = NaN();
    count_passing_matches(
      mcast_matches,
      num_matches, 
      options,
      &num_passing_matches,
      &output_thresh
    );
    fprintf(
      output, 
      "<p>\n"
      "Each of the following %d matches has %s less than %g.\n"
      "<br>\n"
      "The motif matches shown have a position p-value less than %g.\n"
      "<br><b>Click on the arrow</b> (↧) next to the block diagram to view "
      "more information about a sequence.\n"
      "</p>\n",
      num_passing_matches,
      output_thresh_type[options->output_thresh_type],
      output_thresh,
      options->motif_pthresh
    );
  }
  else {
    fprintf(
      output, 
      "<p>\n"
      "MCAST could not evaluate statistics for the match scores. "
      "All the matches found by MCAST are shown below, but no "
      "signficance statistics are provided."
      "</p>\n"
    );
  }

  // Color bar for motifs
  fputs("<div style=\"text-align:left\">\n" ,output);
  int i;
  for (i = 0; i < num_motifs; ++i) {
    MOTIF_T *motif = motif_at(motifs, i);
    MOTIF_T *rc_motif = NULL;
    char *motif_id = get_motif_id(motif);
    char *rc_motif_id = NULL;
    if (i < (num_motifs - 1)) {
      rc_motif = motif_at(motifs, i + 1);
      rc_motif_id = get_motif_id(rc_motif);
    }
    if (rc_motif_id && strcmp(motif_id, rc_motif_id) == 0) {
      ++i; // Pair of identiical motif ids indicate forward/reverse pair.
    }
    fprintf(
      output,
      "<table style=\"display:inline-block; padding:5px;\">"
      "<tr><td>"
      "<div style=\""
      "width: 20px; height: 20px; background-color:%s;border: 1px solid black\">"
      "</div></td>"
      "<td>Motif %s</td>"
      "</tr></table>\n",
      pick_color(i),
      motif_id
    );
  }
  fputs("</div>\n", output);

  // Table of matches
  fputs(
    "<table id=\"tbl_sequences\" style=\"width:100%; table-layout:fixed;\" border=\"0\">\n"
    "<thead><tr>\n"
    "<th class=\"col_more\"> </th>\n"
    "<th class=\"col_seq\">Sequence</th>\n"
    "<th class=\"col_start\">Start</th>\n"
    "<th class=\"col_stop\">Stop</th>\n"
    "<th class=\"col_score\">Score</th>\n"
    "<th class=\"col_pv\">\n"
    "<i>p</i>-value</th>\n"
    "<th class=\"col_ev\">\n"
    "<i>E</i>-value</th>\n"
    "<th class=\"col_qv\">\n"
    "<i>q</i>-value</th>\n"
    "<th class=\"col_extra\"> </th>\n"
    "<th class=\"col_more\"> </th>\n"
    "</tr></thead>\n"
    "<tbody>\n",
    output
  );
  size_t max_seq_len = get_mcast_matches_max_match_length(mcast_matches, num_matches);
  for (i = 0; i < num_passing_matches; ++i) {
    int cluster_id = get_mcast_match_cluster_id(mcast_matches[i]);
    char *seq_name = get_mcast_match_seq_name(mcast_matches[i]);
    size_t start = get_mcast_match_start(mcast_matches[i]);
    size_t stop = get_mcast_match_stop(mcast_matches[i]);
    double score = get_mcast_match_score(mcast_matches[i]);
    double pvalue = get_mcast_match_pvalue(mcast_matches[i]);
    double evalue = get_mcast_match_evalue(mcast_matches[i]);
    double qvalue = get_mcast_match_qvalue(mcast_matches[i]);
    double seq_fract = ((double) (stop - start + 1) / max_seq_len) * 100.0;
    fprintf(
      output,
      "<tr>\n"
      "<td class=\"col_more\"> </td>\n"
      "<td class=\"col_seq\">%s</td>\n"
      "<td class=\"col_start\">%ld</td>\n"
      "<td class=\"col_stop\">%ld</td>\n"
      "<td class=\"col_score\">%g</td>\n"
      "<td class=\"col_pv\">%.3g</td>\n"
      "<td class=\"col_ev\">%.3g</td>\n"
      "<td class=\"col_qv\">%.3g</td>\n"
      "</tr>\n",
      seq_name,
      start,
      stop,
      score,
      pvalue,
      evalue,
      qvalue
    );
    fprintf(
      output,
      "<tr class=\"col_bd\" id=\"cluster-%d_blocks\">\n"
      "<td class=\"col_more\">"
      "<a href=\"javascript:show_more('cluster-%d')\" "
      "class=\"more_arrow\" title=\"Toggle additional information\">↧</a></td>\n"
      "<td class=\"block_td\" colspan=\"8\">"
      "<div class=\"block_container\" id=\"cluster-%d_block_container\">\n"
      "<div class=\"block_plus_sym\">+</div>\n"
      "<div class=\"block_minus_sym\">-</div>\n"
      "<div class=\"block_rule\" style=\"width:%g%%\"></div>\n",
      cluster_id,
      cluster_id,
      cluster_id,
      seq_fract
    );
    MOTIF_HIT_T *motif_hit = retrieve_next_object(mcast_matches[i]->motif_hits);
    while(motif_hit != NULL) {
      fprintf(
        output,
        "<div class=\"block_motif\" style=\"left:%g%%; top:%d; width:%g%%;"
        " height:%g; background-color: %s; border: 1px solid black;\""
        " title=\"%s pvalue: %g starts: %ld ends: %ld\"></div>\n",
        100.0 * ((double) (motif_hit->start - start)) / max_seq_len,
        0,
        100.0 * ((double) (motif_hit->stop - motif_hit->start)) / max_seq_len,
        12.0,
        pick_color(motif_hit->motif_index),
        motif_hit->motif_id,
        motif_hit->pvalue,
        motif_hit->start,
        motif_hit->stop
      );
      motif_hit = retrieve_next_object(mcast_matches[i]->motif_hits);
    }
  }
  fputs(
    "</div></td>\n"
    "</tr>\n"
    "</tbody>\n"
    "</table>\n"
    "</div>\n"
    "</div>\n",
    output
  );
};

void mcast_print_mcast_match_data(
  FILE *output,
  MCAST_MATCH_T **mcast_matches,
  int num_matches,
  MHMMSCAN_OPTIONS_T *options
) {
  int i = 0;
  for (i = 0; i < num_matches; ++i) {
    int cluster_id = mcast_matches[i]->cluster_id;
    size_t start = mcast_matches[i]->start;
    size_t stop = mcast_matches[i]->stop;
    char *sequence = mcast_matches[i]->sequence;
    double pvalue = mcast_matches[i]->pvalue;
    double evalue = mcast_matches[i]->evalue;
    double qvalue = mcast_matches[i]->qvalue;
    if (pvalue <= options->p_thresh 
        && qvalue <= options->q_thresh
        && evalue <= options->e_thresh
    ) {
      fprintf(
        output,
        "<input type=\"hidden\" id=\"cluster-%d_len\" value=\"%ld\"/>\n"
        "<input type=\"hidden\" id=\"cluster-%d_combined_pvalue\" value=\"%3g\"/>\n"
        "<input type=\"hidden\" id=\"cluster-%d_type\" value=\"nucleotide\"/>\n"
        "<input type=\"hidden\" id=\"cluster-%d_segs\" value=\"%ld %s\"/>\n"
        "<input type=\"hidden\" id=\"cluster-%d_hits\" value=\"\n",
        cluster_id,
        stop - start + 1,
        cluster_id,
        pvalue,
        cluster_id,
        cluster_id,
        start,
        sequence,
        cluster_id
      );
      MOTIF_HIT_T *motif_hit = retrieve_next_object(mcast_matches[i]->motif_hits);
      while(motif_hit != NULL) {
        fprintf(
          output,
          "%ld\t%s\t%c\t%3.1g\t%s\n",
          motif_hit->start,
          motif_hit->motif_id,
          motif_hit->strand,
          motif_hit->pvalue,
          motif_hit->seq
        );
        motif_hit = retrieve_next_object(mcast_matches[i]->motif_hits);
      }
      fputs("\"/>\n", output);
    }
  }
}

/**********************************************************************
  mcast_print_results_as_html

  Print an array of mcast_matches as HTML.
**********************************************************************/
void mcast_print_results_as_html(
  MCAST_MATCH_T **mcast_matches,
  int num_matches,
  BOOLEAN stats_available,
  ARRAY_T *bgfreqs,
  MOTIF_T *motifs,
  int num_motifs,
  int num_seqs, 
  long num_residues,
  double duration,
  MHMMSCAN_OPTIONS_T *options
) {
  assert(options != NULL);
  assert(mcast_matches != NULL);

  FILE *mcast_file = fopen(options->html_path, "w");
  if (!mcast_file) {
    die("Couldn't open file %s for output.\n", options->html_path);
  }
  const int MAX_TAG_SIZE = 1000;
  int html_string_size = strlen(mcast_html_string);
  int i = 0;
  for (i = 0; i < html_string_size; ++i) {
    if (mcast_html_string[i] != '@') {
      fputc(mcast_html_string[i], mcast_file);
    }
    else {
      char buffer[MAX_TAG_SIZE];
      ++i;
      int j = 0;
      while (mcast_html_string[i] != '@') {
        buffer[j] = mcast_html_string[i];
        ++j;
        ++i;
      }
      buffer[j] = '\0';
      if (strcmp("sequence_filename", buffer) == 0) {
        mcast_print_sequence_filename(mcast_file, options);
      }
      else if (strcmp("num_sequences", buffer) == 0) {
        mcast_print_num_sequences(mcast_file, num_seqs);
      }
      else if (strcmp("num_residues", buffer) == 0) {
        mcast_print_num_residues(mcast_file, num_residues);
      }
      else if (strcmp("max_match_length", buffer) == 0) {
        mcast_print_max_match_length(mcast_file, mcast_matches, num_matches);
      }
      else if (strcmp("longest_seq_name", buffer) == 0) {
        mcast_print_longest_seq_name(mcast_file, mcast_matches, num_matches);
      }
      else if (strcmp("pattern_filename", buffer) == 0) {
        mcast_print_pattern_filename(mcast_file, options);
      }
      else if (strcmp("motif_array", buffer) == 0) {
        mcast_print_motif_array(mcast_file, motifs, num_motifs);
      }
      else if (strcmp("motif_list", buffer) == 0) {
        mcast_print_motif_list(mcast_file, motifs, num_motifs);
      }
      else if (strcmp("results", buffer) == 0) {
        mcast_print_matches(
          mcast_file,
          mcast_matches,
          num_matches,
          stats_available,
          motifs,
          num_motifs,
          options
        );
      }
      else if (strcmp("version", buffer) == 0) {
        mcast_print_version(mcast_file);
      }
      else if (strcmp("release_date", buffer) == 0) {
        mcast_print_release_date(mcast_file);
      }
      else if (strcmp("command_line", buffer) == 0) {
        mcast_print_command_line(mcast_file, options);
      }
      else if (strcmp("bg_source", buffer) == 0) {
        mcast_print_bg_filename(mcast_file, options);
      }
      else if (strcmp("background_freqs", buffer) == 0) {
        mcast_print_bg_freqs(mcast_file, bgfreqs, options);
      }
      else if (strcmp("mcast_match_data", buffer) == 0) {
        mcast_print_mcast_match_data(
          mcast_file, 
          mcast_matches,
          num_matches,
          options
        );
      }
      else if (strcmp("computation_time", buffer) == 0) {
        mcast_print_duration(mcast_file, duration);
      }
    }
  }
}

/**********************************************************************
  mcast_print_results_as_text

  Print a heap of mcast_matches as plain text.
**********************************************************************/
void mcast_print_results_as_text(
  BOOLEAN_T stats_available,
  int num_matches,
  MCAST_MATCH_T **mcast_matches,
  MHMMSCAN_OPTIONS_T *options
) {

  assert(options != NULL);
  assert(mcast_matches != NULL);

  FILE *out = NULL;
  if (options->text_only) {
    out = stdout;
  }
  else {
    out = fopen(options->text_path, "w");
    if (!out) {
      die("Couldn't open file %s for output.\n", options->text_path);
    }
  }

  // Print header line
  fprintf(
    out,
    "#pattern name\t"
    "sequence name\t"
    "start\t"
    "stop\t"
    "score\t"
    "p-value\t"
    "E-value\t"
    "q-value\t"
    "matched sequence\n"
  );

  int i = 0;
  for (i = 0; i < num_matches; ++i) {
    MCAST_MATCH_T *match = mcast_matches[i];
    double evalue = get_mcast_match_evalue(match);
    double pvalue = get_mcast_match_pvalue(match);
    double qvalue = get_mcast_match_qvalue(match);
    if (!stats_available
       || (
         evalue <= options->e_thresh 
         && pvalue <= options->p_thresh 
         && qvalue <= options->q_thresh)
     ) {
      int id = get_mcast_match_cluster_id(match);
      char *seq_name = get_mcast_match_seq_name(match);
      char *seq = get_mcast_match_sequence(match);
      double score = get_mcast_match_score(match);
      size_t start = get_mcast_match_start(match);
      size_t stop = get_mcast_match_stop(match);
      fprintf(
        out,
        "cluster-%d\t%s\t%zd\t%zd\t%g\t%.5g\t%.5g\t%.5g\t%s\n",
        id,
        seq_name,
        start,
        stop,
        score,
        pvalue,
        evalue,
        qvalue,
        seq+1
      );
    }
  }

  if (out != stdout) {
    fclose(out);
  }
}

