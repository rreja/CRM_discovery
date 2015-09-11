#define _GNU_SOURCE
#include <assert.h>
#include <errno.h>
#include <string.h>
#include <sys/wait.h>

#include "config.h"
#include "dir.h"
#include "eps2png.h"
#include "spamo-output.h"
#include "string-builder.h"
#include "utils.h"
#include "xml-out.h"
#include "seq-reader-from-fasta.h"


#define SPAMO_EPS_TEMPLATE_FILE "spamo_template.eps"

const char* xml_indicator = "<?xml version='1.0' encoding='UTF-8' standalone='yes'?>\n";

const char* spamo_dtd = 
"<!DOCTYPE spamo[\n"
"<!ELEMENT spamo (model, files, primary_motif, run_time)>\n"
"<!ATTLIST spamo version CDATA #REQUIRED release CDATA #REQUIRED>\n"
"<!ELEMENT model (command_line, seed, margin, bin_size,\n"
"  bin_pvalue_calc_range, bin_pvalue_cutoff, seq_max_shared_fract,\n"
"  seq_min_hit_score, redundant_overlap, redundant_joint,\n"
"  motif_pseudocount, motif_trim, bin_max, host, when)>\n"
"<!ELEMENT command_line (#PCDATA)>\n"
"<!ELEMENT seed (#PCDATA)>\n"
"<!ELEMENT margin (#PCDATA)>\n"
"<!ELEMENT bin_size (#PCDATA)>\n"
"<!ELEMENT bin_pvalue_calc_range (#PCDATA)>\n"
"<!ELEMENT bin_pvalue_cutoff (#PCDATA)>\n"
"<!ELEMENT seq_max_shared_fract (#PCDATA)>\n"
"<!ELEMENT seq_min_hit_score (#PCDATA)>\n"
"<!ELEMENT redundant_overlap (#PCDATA)>\n"
"<!ELEMENT redundant_joint (#PCDATA)>\n"
"<!ELEMENT motif_pseudocount (#PCDATA)>\n"
"<!ELEMENT motif_trim (#PCDATA)>\n"
"<!ELEMENT bin_max (#PCDATA)>\n"
"<!ELEMENT host (#PCDATA)>\n"
"<!ELEMENT when (#PCDATA)>\n"
"<!ELEMENT files (sequence_db, motif_db*)>\n"
"<!ELEMENT sequence_db EMPTY>\n"
"<!ATTLIST sequence_db name CDATA #REQUIRED loaded CDATA #REQUIRED\n"
"  excluded_too_short CDATA #REQUIRED excluded_no_match CDATA #REQUIRED\n"
"  excluded_similar CDATA #REQUIRED last_modified CDATA #REQUIRED\n"
"  source CDATA #REQUIRED>\n"
"<!ELEMENT motif_db EMPTY>\n"
"<!ATTLIST motif_db id ID #REQUIRED name CDATA #REQUIRED\n"
"  loaded CDATA #REQUIRED excluded CDATA #REQUIRED\n"
"  last_modified CDATA #REQUIRED source CDATA #REQUIRED\n"
"  cisml CDATA #IMPLIED>\n"
"<!ELEMENT primary_motif (motif, secondary_motif*)>\n"
"<!ELEMENT secondary_motif (spacing*, motif, histogram, redundant?)>\n"
"<!ATTLIST secondary_motif evalue CDATA #REQUIRED>\n"
"<!ELEMENT spacing EMPTY>\n"
"<!ATTLIST spacing strand (same|opposite) #REQUIRED\n"
"  side (left|right) #REQUIRED bin CDATA #REQUIRED num CDATA #REQUIRED\n"
"  pvalue CDATA #REQUIRED>\n"
"<!ELEMENT histogram (same_strand, opposite_strand)>\n"
"<!ATTLIST histogram total CDATA #REQUIRED max CDATA #REQUIRED>\n"
"<!ELEMENT same_strand (left_side, right_side)>\n"
"<!ELEMENT opposite_strand (left_side, right_side)>\n"
"<!ELEMENT left_side (bin*)>\n"
"<!ELEMENT right_side (bin*)>\n"
"<!ELEMENT bin EMPTY>\n"
"<!-- i: index, n: number in bin, p: pvalue (only supplied when bin\n"
"  is tested) -->\n"
"<!ATTLIST bin i CDATA #REQUIRED n CDATA #REQUIRED p CDATA #IMPLIED>\n"
"<!ELEMENT redundant (secondary_motif*)>\n"
"<!-- motif contains the probability of each of the nucleotide bases at each\n"
"  position; i starts at 1; A, C, G and T are probabilities that sum to 1 -->\n"
"<!ELEMENT motif (pos*)>\n"
"<!ATTLIST motif db IDREF #REQUIRED name CDATA #REQUIRED alt CDATA #IMPLIED \n"
"  length CDATA #REQUIRED nsites CDATA #IMPLIED evalue CDATA #IMPLIED\n"
"  ltrim CDATA #IMPLIED rtrim CDATA #IMPLIED url CDATA #IMPLIED>\n"
"<!ELEMENT pos EMPTY>\n"
"<!ATTLIST pos i CDATA #REQUIRED A CDATA #REQUIRED C CDATA #REQUIRED\n"
"  G CDATA #REQUIRED T CDATA #REQUIRED>\n"
"<!-- run time is measured in real time and cpu time -->\n"
"<!ELEMENT run_time EMPTY>\n"
"<!ATTLIST run_time cpu CDATA #REQUIRED real CDATA #REQUIRED>\n"
"]>\n";

#define INT_MAX_WIDTH 15


/**************************************************************************
 * Outputs indent multiples of tab
 * Assumes indent is non-negative
 **************************************************************************/
void output_indent(FILE *file, char *tab, int indent) {
  assert(indent >= 0);
  while (indent--) fputs(tab, file);
}

/**************************************************************************
 * Outputs xml for a fasta database
 **************************************************************************/
void output_sequence_database(FILE *xml_output, SEQUENCE_DB_T *db, char *tab, int indent, char **buffer, int *buffer_len) {
  output_indent(xml_output, tab, indent);
  fputs("<sequence_db ", xml_output);
  
  fprintf(xml_output, "name=\"%s\" loaded=\"%d\" excluded_too_short=\"%d\"\n",
      replace_xml_chars2(db->name, buffer, buffer_len, 0, TRUE), db->loaded, db->excluded_tooshort);

  output_indent(xml_output, tab, indent+2);
  fprintf(xml_output, "excluded_no_match=\"%d\" excluded_similar=\"%d\" last_modified=\"%s\"\n", 
      db->excluded_nomatch, db->excluded_similar, strtok(ctime(&(db->last_mod)), "\n"));

  output_indent(xml_output, tab, indent+2);
  fprintf(xml_output, "source=\"%s\"", replace_xml_chars2(db->source, buffer, buffer_len, 0, TRUE));
  fputs("/>\n", xml_output);
}

/**************************************************************************
 * Outputs xml for a motif database
 **************************************************************************/
void output_motif_database(FILE *xml_output, MOTIF_DB_T *db, char *tab, int indent, char **buffer, int *buffer_len) {
  output_indent(xml_output, tab, indent);
  fputs("<motif_db ", xml_output);
  if (db->id == 0) {
    fputs("id=\"primary_file\"", xml_output);
  } else {
    fprintf(xml_output, "id=\"db%d\"", db->id);
  }
  fprintf(xml_output, " name=\"%s\" loaded=\"%d\" excluded=\"%d\" last_modified=\"%s\"\n",
      replace_xml_chars2(db->name, buffer, buffer_len, 0, TRUE), db->loaded, db->excluded, strtok(ctime(&(db->last_mod)), "\n"));

  output_indent(xml_output, tab, indent+2);
  fprintf(xml_output, "source=\"%s\"", replace_xml_chars2(db->source, buffer, buffer_len, 0, TRUE));
  if (db->cisml) {
    fputs("\n", xml_output);
    output_indent(xml_output, tab, indent+2);
    fprintf(xml_output, "cisml=\"%s\"", replace_xml_chars2(db->cisml, buffer, buffer_len, 0, TRUE));
  }
  fputs("/>\n", xml_output);
}

/**************************************************************************
 * Outputs xml for a motif
 **************************************************************************/
void output_motif(FILE *xml_output, int db_id, MOTIF_T *motif, char *tab, int indent, char **buffer, int *buffer_len) {
  int i, len;
  char *name, *alt, *url;
  MATRIX_T *freqs;
  double A, C, G, T;

  name = get_motif_id(motif);
  alt = get_motif_id2(motif);
  len = get_motif_length(motif);
  url = get_motif_url(motif);
  freqs = get_motif_freqs(motif);

  output_indent(xml_output, tab, indent);
  fputs("<motif ", xml_output);
  if (db_id == 0) {
    fputs("db=\"primary_file\"", xml_output);
  } else {
    fprintf(xml_output, "db=\"db%d\"", db_id);
  }
  fprintf(xml_output, " name=\"%s\" ", replace_xml_chars2(name, buffer, buffer_len, 0, TRUE));
  if (alt && alt[0] != '\0') fprintf(xml_output, "alt=\"%s\" ", replace_xml_chars2(alt, buffer, buffer_len, 0,  TRUE));
  fprintf(xml_output, "length=\"%d\" nsites=\"%g\" evalue=\"%g\" ", len, get_motif_nsites(motif), get_motif_evalue(motif));
  fprintf(xml_output, "ltrim=\"%d\" rtrim=\"%d\" ", get_motif_trim_left(motif), get_motif_trim_right(motif));
  if (url && url[0] != '\0') fprintf(xml_output, "url=\"%s\"", replace_xml_chars2(url, buffer, buffer_len, 0,  TRUE));
  fprintf(xml_output, ">\n");
  for (i = 0; i < len; ++i) {
    A = get_matrix_cell(i, 0, freqs);
    C = get_matrix_cell(i, 1, freqs);
    G = get_matrix_cell(i, 2, freqs);
    T = get_matrix_cell(i, 3, freqs);
    output_indent(xml_output, tab, indent + 1);
    fprintf(xml_output, "<pos i=\"%d\" A=\"%g\" C=\"%g\" G=\"%g\" T=\"%g\"/>\n", (i+1), A, C, G, T);
  }
  output_indent(xml_output, tab, indent);
  fprintf(xml_output, "</motif>\n");
}

/**************************************************************************
 * Outputs xml for the histogram spacings for one quadrant.
 **************************************************************************/
static void output_bins(FILE *xml_output, SPACING_T *spacings, int bin_count, char *tab, int indent) {
  int i;
  for (i = 0; i < bin_count; ++i) {
    output_indent(xml_output, tab, indent);
    if (spacings->pvalues[i] <= 1 && spacings->bins[i] > 0) { 
      fprintf(xml_output, "<bin i=\"%d\" n=\"%d\" p=\"%.2g\"/>\n", i+1, spacings->bins[i], spacings->pvalues[i]);
    } else {
      fprintf(xml_output, "<bin i=\"%d\" n=\"%d\"/>\n", i+1, spacings->bins[i]);
    }
  }
}

/**************************************************************************
 * Outputs xml for the histogram of the spacings
 **************************************************************************/
static void output_histogram(FILE *xml_output, SECONDARY_MOTIF_T *smotif, int margin, int bin, char *tab, int indent) {
  int max, possible, bin_count;
  possible = margin - get_motif_trimmed_length(smotif->motif) + 1;
  bin_count = (possible / bin) + (possible % bin ? 1 : 0);
  max = smotif->max_in_one_bin;
  output_indent(xml_output, tab, indent);
  fprintf(xml_output, "<histogram total=\"%d\" max=\"%d\">\n", smotif->total_spacings, max);
  output_indent(xml_output, tab, indent + 1);
  fprintf(xml_output, "<same_strand>\n");
  output_indent(xml_output, tab, indent + 2);
  fprintf(xml_output, "<left_side>\n");
  output_bins(xml_output, smotif->spacings+(SAME | LEFT), bin_count, tab, indent + 3);
  output_indent(xml_output, tab, indent + 2);
  fprintf(xml_output, "</left_side>\n");
  output_indent(xml_output, tab, indent + 2);
  fprintf(xml_output, "<right_side>\n");
  output_bins(xml_output, smotif->spacings+(SAME | RIGHT), bin_count, tab, indent + 3);
  output_indent(xml_output, tab, indent + 2);
  fprintf(xml_output, "</right_side>\n");
  output_indent(xml_output, tab, indent + 1);
  fprintf(xml_output, "</same_strand>\n");
  output_indent(xml_output, tab, indent + 1);
  fprintf(xml_output, "<opposite_strand>\n");
  output_indent(xml_output, tab, indent + 2);
  fprintf(xml_output, "<left_side>\n");
  output_bins(xml_output, smotif->spacings+(OPPO | LEFT), bin_count, tab, indent + 3);
  output_indent(xml_output, tab, indent + 2);
  fprintf(xml_output, "</left_side>\n");
  output_indent(xml_output, tab, indent + 2);
  fprintf(xml_output, "<right_side>\n");
  output_bins(xml_output, smotif->spacings+(OPPO | RIGHT), bin_count, tab, indent + 3);
  output_indent(xml_output, tab, indent + 2);
  fprintf(xml_output, "</right_side>\n");
  output_indent(xml_output, tab, indent + 1);
  fprintf(xml_output, "</opposite_strand>\n");
  output_indent(xml_output, tab, indent);
  fprintf(xml_output, "</histogram>\n");
}

/**************************************************************************
 * Outputs xml for the histogram of the spacings
 **************************************************************************/
void output_secondary_motif(FILE *xml_output, SECONDARY_MOTIF_T *smotif, 
    int n_secondary_motifs, LINKLST_T *rmotifs, 
    int margin, int bin, char *tab, int indent, char **buffer, int *buffer_len) {
  LINK_T *node;
  SIGSPACE_T sig;
  int i;
  output_indent(xml_output, tab, indent);
  fprintf(xml_output, "<secondary_motif evalue=\"%.2g\">\n", smotif->min_pvalue * n_secondary_motifs);
  output_indent(xml_output, tab, indent + 1);
  for (i = 0; i < smotif->sig_count; ++i) {
    sig = smotif->sigs[i];
    output_indent(xml_output, tab, indent + 1);
    fprintf(xml_output, "<spacing strand=\"%s\" side=\"%s\" bin=\"%d\" num=\"%d\" pvalue=\"%.2g\"/>\n", 
        (STRAND(sig.quad) == SAME ? "same" : "opposite" ),
        (SIDE(sig.quad) == LEFT ? "left" : "right" ), 
        sig.bin + 1, smotif->spacings[sig.quad].bins[sig.bin], sig.pvalue);
  }
  output_motif(xml_output, smotif->db->id, smotif->motif, tab, indent + 1, buffer, buffer_len);
  output_histogram(xml_output, smotif, margin, bin, tab, indent + 1);
  if (rmotifs != NULL && linklst_size(rmotifs) > 0) {
    output_indent(xml_output, tab, indent + 1);
    fputs("<redundant>\n", xml_output);
    for (node = linklst_first(rmotifs); node != NULL; node = linklst_next(node)) {
      output_secondary_motif(xml_output, (SECONDARY_MOTIF_T*)linklst_get(node), n_secondary_motifs,
        NULL, margin, bin, tab, indent + 2, buffer, buffer_len);
    }
    output_indent(xml_output, tab, indent + 1);
    fputs("</redundant>\n", xml_output);
  }
  output_indent(xml_output, tab, indent);
  fprintf(xml_output, "</secondary_motif>\n");
}

/**************************************************************************
 * Outputs an encapsulated postscript histogram representing the spacings 
 * of the secondary motif to the primary.
 **************************************************************************/
static BOOLEAN_T output_eps(FILE *file, int margin, int binsize, int binmax, double sigthresh, FILE *template, SECONDARY_MOTIF_T *smotif) {
  char *title;
  int ch, len, sequences, opts, count, i;
  SPACING_T *spacing;
  //copy the template to the destination file
  rewind(template);
  while (1) {
    if ((ch = fgetc(template)) == EOF) {
      if (ferror(template)) return FALSE;
      break;
    }
    if (fputc(ch, file) == EOF) return FALSE;
  }
  title = "";
  len = get_motif_trimmed_length(smotif->motif);
  opts = margin - len + 1;
  count = (opts / binsize) + (opts % binsize ? 1 : 0);
  sequences = smotif->total_spacings;

  //write out the graph data
  if (
    fprintf(file, 
        "\n\n%% Input Data\n"
        "(%s) title\n"
        "%d margin\n"
        "%d motif-length\n"
        "%d bin-size\n"
        "%d bin-max\n"
        "%d sequences\n"
        "%g threshold\n",
        title, margin, len, binsize, binmax, sequences, sigthresh) < 0
   ) return FALSE;

  spacing = &(smotif->spacings[SAME | LEFT]);
  if (fputs("\n\n% Same Left\n", file) == EOF) return FALSE;
  for (i = 0; i < count; ++i) {
    if (fprintf(file, "%d %g\n", spacing->bins[i], spacing->pvalues[i]) < 0) return FALSE;
  }
  spacing = &(smotif->spacings[SAME | RIGHT]);
  if (fputs("\n\n% Same Right\n", file) == EOF) return FALSE;
  for (i = 0; i < count; ++i) {
    if (fprintf(file, "%d %g\n", spacing->bins[i], spacing->pvalues[i]) < 0) return FALSE;
  }
  spacing = &(smotif->spacings[OPPO | LEFT]);
  if (fputs("\n\n% Opposite Left\n", file) == EOF) return FALSE;
  for (i = 0; i < count; ++i) {
    if (fprintf(file, "%d %g\n", spacing->bins[i], spacing->pvalues[i]) < 0) return FALSE;
  }
  spacing = &(smotif->spacings[OPPO | RIGHT]);
  if (fputs("\n\n% Opposite Right\n", file) == EOF) return FALSE;
  for (i = 0; i < count; ++i) {
    if (fprintf(file, "%d %g\n", spacing->bins[i], spacing->pvalues[i]) < 0) return FALSE;
  }
  if (fputs("\n\ngraph\nshowpage", file) == EOF) return FALSE;
  return TRUE;
}

/**************************************************************************
 * Makes the histogram file name from the details in the motifs and the
 * file extension.
 * Caller is responsible for freeing memory
 **************************************************************************/
static char* make_pattern_file_name(char *prefix, char *ext, MOTIF_T* primary, SECONDARY_MOTIF_T *secondary) {
  const char *fmt = "%s_%s_db%d_%s.%s";
  char dummy[1];
  char *ret;
  int len;
  len = snprintf(dummy, 1, fmt, prefix, get_motif_id(primary), secondary->db->id, get_motif_id(secondary->motif), ext);
  ret = mm_malloc(sizeof(char) * (len+1));
  snprintf(ret, (len+1), fmt, prefix, get_motif_id(primary), secondary->db->id, get_motif_id(secondary->motif), ext);
  return ret;
}


/**************************************************************************
 * Creates a histogram
 **************************************************************************/
static void create_histogram
(int margin, int binsize, int binmax, double sigthresh, char *dir, 
 FILE *template, MOTIF_T* primary, SECONDARY_MOTIF_T *secondary, 
 BOOLEAN_T make_eps, BOOLEAN_T make_png)
{
  STR_T *path;
  FILE *file;
  path = str_create(0);
  // check that it's not the current directory
  if (strcmp(dir, "") != 0 && strcmp(dir, ".") != 0) { 
    str_append(path, dir, strlen(dir));
    if (str_char(path, -1) != '/') str_append(path, "/", 1);
  }
  str_appendf(path, "hist_%s_db%d_%s", get_motif_id(primary), 
      secondary->db->id, get_motif_id(secondary->motif));
  
  // make the eps
  str_append(path, ".eps", 4);
  file = fopen(str_internal(path), "w");
  output_eps(file, margin, binsize, binmax, sigthresh, template, secondary);
  fclose(file);
  str_truncate(path, -4);

  // make the png
  if (make_png) eps_to_png(str_internal(path));

  // possibly remove the eps
  if (!make_eps) {
    str_append(path, ".eps", 4);
    unlink(str_internal(path));
  }

  // cleanup
  str_destroy(path, FALSE);
}

/**************************************************************************
 * Creates histograms
 **************************************************************************/
void create_histograms
(int margin, int binsize, int binmax, double sigthresh, char *dir, 
 MOTIF_T *primary_motif, LINKLST_T *secondary_motifs, 
 BOOLEAN_T make_eps, BOOLEAN_T make_png, BOOLEAN_T make_for_redundant) 
{
  char *template_path;
  FILE *template;
  LINK_T *node;
  GROUPED_MOTIF_T *gmotif;
  //check that we can convert eps to png (if png wanted)
  if (make_png && !eps_to_png_possible()) {
    fprintf(stderr, "Warning: No usable eps conversion programs. Skiping creation of png histograms.\n");
    make_png = FALSE; //can't create a png without ghostscript or convert
  }
  //check that histograms are actually wanted
  if (!(make_eps || make_png)) return;
  //check that we have access to the template, otherwise we can't create any histograms
  template_path = make_path_to_file(get_meme_etc_dir(), SPAMO_EPS_TEMPLATE_FILE);
  if (!file_exists(template_path)) {
    fprintf(stderr, "Warning: Can not find SpaMo histogram template file \"%s\". Skiping creation of histograms.\n", template_path);
    free(template_path);
    return;
  }

  //open the template
  template = fopen(template_path, "r");
  //loop over the non-redundant significant motifs and make histograms
  for (node = linklst_first(secondary_motifs); node != NULL; node = linklst_next(node)) {
    gmotif = (GROUPED_MOTIF_T*)linklst_get(node);
    create_histogram(margin, binsize, binmax, sigthresh, dir, template, 
        primary_motif, gmotif->best, make_eps, make_png);
    //optionally make the histograms for the redundant motifs as well
    if (make_for_redundant) {
      LINK_T *other;
      for (other = linklst_first(gmotif->others); other != NULL; 
          other = linklst_next(other)) {
        create_histogram(margin, binsize, binmax, sigthresh, dir, template, 
            primary_motif, (SECONDARY_MOTIF_T*)linklst_get(other), 
            make_eps, make_png);
      }
    }
  }
  fclose(template);
  free(template_path);
}


/**************************************************************************
 * Dump sequence matches sorted by the name of the sequence.
 *
 * Outputs Columns:
 *   1) Trimmed lowercase sequence with uppercase matches.
 *   2) Position of the secondary match within the whole sequence.
 *   3) Sequence fragment that the primary matched.
 *   4) Strand of the primary match (+|-)
 *   5) Sequence fragment that the secondary matched.
 *   6) Strand of the secondary match (+|-)
 *   7) Is the primary match on the same strand as the secondary (s|o)
 *   8) Is the secondary match downstream or upstream (d|u)
 *   9) The gap between the primary and secondary matches
 *  10) The name of the sequence
 *  11) The p-value of the bin containing the match (adjusted for # of bins)
 *  ---if the FASTA input file sequence names are in Genome Browser format:
 *  12-14) Position of primary match in BED coordinates
 *  15) Position of primary match in Genome Browser coordinates
 *  16-18) Position of secondary match in BED coordinates
 *  19) Position of secondary match in Genome Browser coordinates
 *
 * If you wish to sort based on the gap column:
 * Sort individual output:
 *  sort -n -k 9,9 -o seqs_primary_secondary.txt seqs_primary_secondary.txt
 * Or sort all outputs:
 *  for f in seqs_*.txt; do sort -n -k 9,9 -o $f $f; done
 * Or to get just locations of primary motif in BED coordinates
 * where the secondary is on the opposite strand, upstream with a gap of 118bp:
 *   awk '$7=="o" && $8=="u" && $9==118 {print $12"\t"$13"\t"$14;}' seqs_primary_secondary.txt 
 *
 **************************************************************************/
void dump_sequence_matches(FILE *out, int margin, int bin, double sigthresh, BOOLEAN_T sig_only, RBTREE_T *sequences, MOTIF_T *primary_motif, SECONDARY_MOTIF_T *secondary_motif, int *matches) {
  RBNODE_T *node;
  SEQUENCE_T *sequence;
  int seqlen, i, j, start, end, secondary, secondary_pos, primary_len, secondary_len, distance;
  BOOLEAN_T primary_rc, secondary_rc, downstream; 
  char *buffer, *seq, *primary_match, *secondary_match;
  // allocate a buffer for copying the trimmed sequence into and modify it
  seqlen = margin * 2 + get_motif_trimmed_length(primary_motif);
  buffer = (char*)mm_malloc(sizeof(char) * (seqlen + 1));
  // get the lengths of the motifs
  primary_len = get_motif_trimmed_length(primary_motif);
  secondary_len = get_motif_trimmed_length(secondary_motif->motif); 
  // allocate some strings for storing the matches
  primary_match = (char*)mm_malloc(sizeof(char) * (primary_len + 1));
  secondary_match = (char*)mm_malloc(sizeof(char) * (secondary_len + 1));
  // add null byte at the end of the match strings
  primary_match[primary_len] = '\0';
  secondary_match[secondary_len] = '\0';

  // iterate over all the sequences
  for (node = rbtree_first(sequences); node != NULL; node = rbtree_next(node)) {
    sequence = (SEQUENCE_T*)rbtree_value(node);
    primary_rc = sequence->primary_match < 0;
    secondary = matches[sequence->index];
    if (secondary == 0) continue;
    secondary_rc = secondary < 0;
    secondary_pos = abs(secondary);

    // calculate the distance
    if (secondary_pos <= margin) {
      distance = margin - secondary_pos - secondary_len + 1;
      downstream = primary_rc;
    } else {
      distance = secondary_pos - margin - primary_len - 1;
      downstream = !primary_rc;
    }

    // copy the trimmed sequence
    seq = sequence->data;
    for (i = 0; i < seqlen; ++i) {
      buffer[i] = tolower(seq[i]);
    }
    buffer[seqlen] = '\0';

    // uppercase primary
    start = margin;
    end = margin + primary_len;
    for (i = start, j = 0; i < end; ++i, ++j) {
      buffer[i] = toupper(buffer[i]);
      primary_match[j] = buffer[i];
    }

    // uppercase secondary
    // note orign was one, subtract 1 to make origin zero as required for arrays
    start = secondary_pos -1;
    end = start + secondary_len;
    for (i = start, j = 0; i < end; ++i, ++j) {
      buffer[i] = toupper(buffer[i]);
      secondary_match[j] = buffer[i];
    }

    // get the p-value of the seconndary match
    SPACING_T *spacings;
    if (secondary_rc == primary_rc) {
      spacings = downstream ? secondary_motif->spacings+(SAME+RIGHT) : secondary_motif->spacings+(SAME+LEFT); 
    } else {
      spacings = downstream ? secondary_motif->spacings+(OPPO+RIGHT) : secondary_motif->spacings+(OPPO+LEFT); 
    }
    double p_value = spacings->pvalues[distance/bin];

    // skip match if not significant and only reporting significant matches
    if (sig_only && (p_value > sigthresh)) continue;

    // output line to file
    fprintf(out, "%s    %3d    %s    %s    %s    %s    %s    %s    %3d    %s    %.1e", 
        buffer, 
        secondary_pos, 
        primary_match, 
        (primary_rc ? "-" : "+"), 
        secondary_match, 
        (secondary_rc ? "-" : "+"), 
        (secondary_rc == primary_rc ? "s" : "o"),
        (downstream ? "d" : "u"), 
        distance, 
        sequence->name,
        p_value
    );

    // Parse the sequence name to see if we can get genomic coordinates
    // and print additional columns with primary and secondary matches
    // in both BED and Genome Browser coordinates.
    char *chr_name;
    size_t chr_name_len;
    int start_pos, end_pos;
    if (parse_genomic_coordinates_helper(
        sequence->name,
        &chr_name,
	&chr_name_len,
	&start_pos,
        &end_pos))
     {
       // Get the start and end of the primary match in 
       // 0-relative, half-open genomic coordinates.
       int p_start = start_pos + abs(sequence->primary_match) - 1;
       int p_end = p_start + primary_len;
       // Get the start and end of the secondary match in 
       // 0-relative, half-open genomic coordinates.
       int s_start, s_end;
       if ( (!primary_rc && downstream) || (primary_rc && !downstream) ) {
         s_start = p_end + distance;
         s_end = s_start + secondary_len;
       } else {
         s_end = p_start - distance;
         s_start = s_end - secondary_len;
       }
       fprintf(out, "    %s    %d    %d    %s:%d-%d", 
         chr_name, p_start, p_end, chr_name, p_start+1, p_end);
       fprintf(out, "    %s    %d    %d    %s:%d-%d\n", 
         chr_name, s_start, s_end, chr_name, s_start+1, s_end);
     } else {
       fprintf(out, "\n");
     }
  }

  free(buffer);
  free(primary_match);
  free(secondary_match);
}

/**************************************************************************
 * Create an output file and dump the sequence matches to file.
 **************************************************************************/
void output_sequence_matches(char *dir, int margin, int bin, double sigthresh, BOOLEAN_T sig_only, RBTREE_T *sequences, MOTIF_T *primary_motif, SECONDARY_MOTIF_T *secondary_motif, int *matches) {
  FILE *out;
  int file_name_len;
  char *file_path, *file_name;
  file_name = make_pattern_file_name("seqs", "txt", primary_motif, secondary_motif);
  file_path = make_path_to_file(dir, file_name);
  out = fopen(file_path, "w");
  dump_sequence_matches(out, margin, bin, sigthresh, sig_only, sequences, primary_motif, secondary_motif, matches);
  fclose(out);
  free(file_path);
  free(file_name);
}

