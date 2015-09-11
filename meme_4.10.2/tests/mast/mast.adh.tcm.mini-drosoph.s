<?xml version='1.0' encoding='UTF-8' standalone='yes'?>
<!DOCTYPE mast[
<!ELEMENT mast (model, alphabet, motifs, sequences, runtime)>
<!ATTLIST mast version CDATA #REQUIRED release CDATA #REQUIRED>
<!ELEMENT model (command_line, max_correlation, remove_correlated, strand_handling, translate_dna, max_seq_evalue,
    adj_hit_pvalue, max_hit_pvalue, max_weak_pvalue, host, when)>
<!ELEMENT command_line (#PCDATA)>
<!ELEMENT max_correlation (#PCDATA)>
<!ELEMENT remove_correlated EMPTY>
<!ATTLIST remove_correlated value (y|n) #REQUIRED>
<!ELEMENT strand_handling EMPTY>
<!ATTLIST strand_handling value (combine|separate|norc|protein) #REQUIRED>
<!ELEMENT translate_dna EMPTY>
<!ATTLIST translate_dna value (y|n) #REQUIRED>
<!ELEMENT max_seq_evalue (#PCDATA)>
<!ELEMENT adj_hit_pvalue EMPTY>
<!ATTLIST adj_hit_pvalue value (y|n) #REQUIRED>
<!ELEMENT max_hit_pvalue (#PCDATA)>
<!ELEMENT max_weak_pvalue (#PCDATA)>
<!ELEMENT host (#PCDATA)>
<!ELEMENT when (#PCDATA)>
<!ELEMENT alphabet (letter*)>
<!ATTLIST alphabet type (amino-acid|nucleotide) #REQUIRED bg_source (preset|file|sequence_composition) #REQUIRED bg_file CDATA #IMPLIED>
<!ELEMENT letter EMPTY>
<!ATTLIST letter symbol CDATA #REQUIRED ambig (y|n) "n" bg_value CDATA #IMPLIED>
<!ELEMENT motifs (motif*,correlation*,nos*)>
<!ATTLIST motifs source CDATA #REQUIRED name CDATA #REQUIRED last_mod_date CDATA #REQUIRED>
<!ELEMENT motif EMPTY>
<!-- num is simply the loading order of the motif, it's superfluous but makes things easier for XSLT -->
<!ATTLIST motif id ID #REQUIRED num CDATA #REQUIRED name CDATA #REQUIRED width CDATA #REQUIRED
   best_f CDATA #REQUIRED best_r CDATA #IMPLIED bad (y|n) "n">
<!-- for n > 1 motifs there should be (n * (n - 1)) / 2 correlations, obviously there are none for only 1 motif -->
<!ELEMENT correlation EMPTY>
<!ATTLIST correlation motif_a IDREF #REQUIRED motif_b IDREF #REQUIRED value CDATA #REQUIRED>
<!-- nos: Nominal Order and Spacing diagram, a rarely used feature where mast can adjust pvalues for an expected motif spacing -->
<!ELEMENT nos (expect*)>
<!-- length is in the same unit as the motifs, which is not always the same unit as the sequence -->
<!ATTLIST nos length CDATA #REQUIRED>
<!-- the expect tags are expected to be ordered by pos ascending -->
<!ELEMENT expect EMPTY>
<!ATTLIST expect pos CDATA #REQUIRED gap CDATA #REQUIRED motif IDREF #REQUIRED>
<!ELEMENT sequences (database*, sequence*)>
<!-- the database tags are expected to be ordered in file specification order -->
<!ELEMENT database EMPTY>
<!ATTLIST database id ID #REQUIRED num CDATA #REQUIRED source CDATA #REQUIRED name CDATA #REQUIRED last_mod_date CDATA #REQUIRED 
    seq_count CDATA #REQUIRED residue_count CDATA #REQUIRED type (amino-acid|nucleotide) #REQUIRED link CDATA #IMPLIED>
<!-- the sequence tags are expected to be ordered by best combined p-value (of contained score tags) ascending -->
<!ELEMENT sequence (score*,seg*)>
<!ATTLIST sequence id ID #REQUIRED db IDREF #REQUIRED num CDATA #REQUIRED name CDATA #REQUIRED comment CDATA "" length CDATA #REQUIRED>
<!ELEMENT score EMPTY>
<!-- frame is the starting offset for translation of dna sequences which gives the lowest pvalues for the provided protein motifs -->
<!ATTLIST score strand (both|forward|reverse) #REQUIRED frame (a|b|c) #IMPLIED combined_pvalue CDATA #REQUIRED evalue CDATA #REQUIRED>
<!-- within each sequence the seg tags are expected to be ordered by start ascending -->
<!ELEMENT seg (data,hit*)>
<!ATTLIST seg start CDATA #REQUIRED>
<!ELEMENT data (#PCDATA)>
<!-- within each seg the hit tags are expected to be ordered by pos ascending and then forward strand first -->
<!ELEMENT hit EMPTY>
<!-- gap, while superfluous, makes creating motif diagrams for the text version much easier when using XSLT -->
<!ATTLIST hit pos CDATA #REQUIRED gap CDATA #REQUIRED motif IDREF #REQUIRED pvalue CDATA #REQUIRED strand (forward|reverse) "forward" 
    match CDATA #REQUIRED translation CDATA #IMPLIED>
<!ELEMENT runtime EMPTY>
<!ATTLIST runtime cycles CDATA #REQUIRED seconds CDATA #REQUIRED>
]>
<mast version="4.8.1" release="Mon Jun 25 12:20:25 EST 2012">
	<model>
		<command_line>mast /home/james/meme/scripts/../tests/meme/meme.adh.tcm - -oc /home/james/meme/scripts/../tests/results/mast -nohtml -nostatus -df stdin -dna -seqp</command_line>
		<max_correlation>0.60</max_correlation>
		<remove_correlated value="n"/>
		<strand_handling value="combine"/>
		<translate_dna value="y"/>
		<max_seq_evalue>10</max_seq_evalue>
		<adj_hit_pvalue value="y"/>
		<max_hit_pvalue>0.0001</max_hit_pvalue>
		<max_weak_pvalue>0.0001</max_weak_pvalue>
		<host>tbl-squirrel</host>
		<when>Mon Jun 25 12:50:46 2012</when>
	</model>
	<alphabet type="amino-acid" bg_source="preset">
		<letter symbol="A" bg_value="0.070"/>
		<letter symbol="C" bg_value="0.024"/>
		<letter symbol="D" bg_value="0.040"/>
		<letter symbol="E" bg_value="0.052"/>
		<letter symbol="F" bg_value="0.040"/>
		<letter symbol="G" bg_value="0.074"/>
		<letter symbol="H" bg_value="0.029"/>
		<letter symbol="I" bg_value="0.041"/>
		<letter symbol="K" bg_value="0.052"/>
		<letter symbol="L" bg_value="0.096"/>
		<letter symbol="M" bg_value="0.017"/>
		<letter symbol="N" bg_value="0.032"/>
		<letter symbol="P" bg_value="0.065"/>
		<letter symbol="Q" bg_value="0.042"/>
		<letter symbol="R" bg_value="0.067"/>
		<letter symbol="S" bg_value="0.084"/>
		<letter symbol="T" bg_value="0.052"/>
		<letter symbol="V" bg_value="0.059"/>
		<letter symbol="W" bg_value="0.016"/>
		<letter symbol="Y" bg_value="0.022"/>
	</alphabet>
	<motifs source="/home/james/meme/scripts/../tests/meme/meme.adh.tcm" name="/home/james/meme/scripts/../tests/meme/meme.adh.tcm" last_mod_date="Fri Jun 22 17:44:15 2012">
		<motif id="motif_1" num="1" name="1" width="21" best_f="LITGCSSGIGKAIAKHFHKEG"/>
		<motif id="motif_2" num="2" name="2" width="22" best_f="YCASKFAVRGFTRSMAMEYAPY"/>
		<correlation motif_a="motif_1" motif_b="motif_2" value="0.34"/>
	</motifs>
	<sequences>
		<database id="db_1" num="1" source="-" name="stdin" last_mod_date="Mon Jun 25 12:50:46 2012" seq_count="4" residue_count="499297" type="nucleotide"/>
		<sequence id="seq_1_3" db="db_1" num="3" name="gi|7296683|gb|AE003602.1|AE003602" comment="Drosophila melanogaster genomic scaffold 142000013386043 section 3 of 8, complete sequence" length="297266">
			<score strand="both" combined_pvalue="1.57e-16" evalue="6.3e-16" frame="b"/>
			<seg start="101626">
				<data>
CCCTTACATAATTCTGGATCGGATCATGAATTTCGCGGGCAAAGTGGTCCTTATTACGGGAGCAAGCTCCGGAAT
CGGAGCTGCAACCGCCATTAAGTTTGCCAAGTACGGCGCCTGTCTGGCTCTCAATGGACGCAATGTGGAGAACCT
				</data>
				<hit pos="101675" gap="101674" motif="motif_1" pvalue="3.4e-11" strand="forward" match="++++++++++ ++++ +++++" translation="LITGASSGIGAATAIKFAKYG"/>
			</seg>
			<seg start="102901">
				<data>
GCAAGGGCAAGGCTGGTCTGGACTTCTCCGGCAAAGTGGTGCTTATCACGGGCGCAGCCTCCGGGATCGGGGCCG
CCGCGGCGGAGATGTTCTCGAAGCTGGGTGCCTGCCTGGCCCTGGTGGATCGGGAGGAGGAGGGCCTCATATGTG
				</data>
				<hit pos="102942" gap="1204" motif="motif_1" pvalue="1.8e-08" strand="forward" match="++++++++++ + ++ + + +" translation="LITGAASGIGAAAAEMFSKLG"/>
			</seg>
			<seg start="103351">
				<data>
TTCGCGCTTTTCCCAATCTGGTGGCCTACAACATGTCCAAGGCGGCGGTGGACCAGTTTACCCGCTCCCTTGCCC
TGGATCTGGGTCCCCAGGGTGTTCGGGTGAATGCGGTCAATCCGGGTGTGATTCGCACCAATCTGCAAAAGGCGG
				</data>
				<hit pos="103377" gap="372" motif="motif_2" pvalue="2.8e-08" strand="forward" match="+ +++++++ ++++++++++++" translation="YNMSKAAVDQFTRSLALDLGPQ"/>
			</seg>
			<seg start="104176">
				<data>
CACCACGACGCTGCAGCTCGGTGATGATTACGCCGGGATTCACGGAGTTCACACGCACACCCTTGGGAGCTAGCT
CCAGAGCCACGCACCTGGTGAACTGATCCACGGCAGCCTTGGAAACATTGTATGCTAAGACTCCGGGAAAGGAAC
				</data>
				<hit pos="104237" gap="794" motif="motif_2" pvalue="2.7e-07" strand="reverse" match="+++++++ ++++ ++++++  +" translation="KPALELAVCRTFQDVAAKSVNY"/>
			</seg>
			<seg start="104626">
				<data>
CGGTCTCGTTGAGCTTATCCAAATTCCTGCCCACGATGGTGAGCAGGCCTCCCAGTTTAGCCAAGAGCACCGAAG
TACCCGCTCCAATTCCCGAACTGGCTCCGGTCACGATTATAACTTTATCCTTGAAACTGGGCATGTTGAAGTGAT
				</data>
				<hit pos="104675" gap="372" motif="motif_1" pvalue="4.5e-08" strand="reverse" match="+ +++  ++  ++++++++++" translation="GLKALLVSTGAGIGSSAGTVI"/>
			</seg>
			<seg start="105526">
				<data>
GGACTCCTTTGGGGGCCAGTTCCAGGGCTATGCAGGCCGTGAACTGGTCCACCGCTGCCTTAGACACATTGTAGG
				</data>
				<hit pos="105533" gap="795" motif="motif_2" pvalue="6.2e-06" strand="reverse" match="+++++++ + ++ ++++++  +" translation="KPALELAICATFQDVAAKSVNY"/>
			</seg>
			<seg start="105901">
				<data>
CGGCCGCCACAATGTTGTCGGCCGTCTCCTTCAGCTTCTCCTCGTTGCGACCCACGATGACCAGGAGCCCTCCGA
GTTTCGCCAAATGGACGGCGGCACTTGCTCCGATGCCGGAGCTGGCGCCGGTAACAATAATCACTTTGTCCTTGA
				</data>
				<hit pos="105971" gap="372" motif="motif_1" pvalue="1.8e-09" strand="reverse" match="+ ++++ + + ++++++++++" translation="GLKALHVAASAGIGSSAGTVI"/>
			</seg>
		</sequence>
		<sequence id="seq_1_4" db="db_1" num="4" name="gi|7295475|gb|AE003567.1|AE003567" comment="Drosophila melanogaster genomic scaffold 142000013386050 section 54 of 54, complete sequence" length="104808">
			<score strand="both" combined_pvalue="2.33e-07" evalue="9.3e-07" frame="c"/>
			<seg start="87151">
				<data>
CCTGGAGCCAGGCAGTTGACGCGAATGCCCTCGGGCGCCAGATCCTTGGCGGCTGCCTTGGTCAAGCCAATCAGC
GCGGTCTTGCTGACGGAATAGGCTCCCAGTAGCTATGCGATTAAGATAACGGAGATAAGCATTGAACAACTGAAC
				</data>
				<hit pos="87180" gap="87179" motif="motif_2" pvalue="3.0e-08" strand="reverse" match="++++++++++++++++ ++ ++" translation="EPALDKAAAKTLGILATKSVSY"/>
			</seg>
		</sequence>
		<sequence id="seq_1_2" db="db_1" num="2" name="gi|7289065|gb|AE002569.1|AE002569" comment="Drosophila melanogaster genomic scaffold 142000013385354, complete sequence" length="12850">
			<score strand="both" combined_pvalue="1.01e-01" evalue="0.41" frame="a"/>
		</sequence>
		<sequence id="seq_1_1" db="db_1" num="1" name="gi|7289301|gb|AE002567.1|AE002567" comment="Drosophila melanogaster genomic scaffold 142000013385554, complete sequence" length="84373">
			<score strand="both" combined_pvalue="3.71e-01" evalue="1.5" frame="a"/>
		</sequence>
	</sequences>
	<runtime cycles="1448090" seconds="1.448"/>
</mast>
