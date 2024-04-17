#!/bin/bash

# change to use Vertebrata_OrthoDB_11 db instead of vertebrate_aa which was OrthoDB version 10

#### define global vars #####################################################################################################

# these files are soft links to actual files. soft links created when running functional_anno_setup.sh
in_gff=annot.gff
in_faa_file=proteins.faa
in_cds_file=codingseq.fna

gff=$in_gff
protein_file=$in_faa_file
cds_file=$in_cds_file

# interproscan vars
prefix=ipr/interproscan
ipr_tsv=${prefix}.tsv
ipr_anno=${prefix}.anno
ipr_end_file=${prefix}_end

tsv_dir=blasts

nt_tsv=nt_blastn.tsv

unmatched_cds_file=unmatched_by_blastp.codingseq.fna
unmatched_nt_tsv=nt_run_on_unmatched_recs_blastn.tsv

unmatched_protein_file=unmatched_by_tsvs.faa
umatched_protein_tsv=blastp_on_unmatched.tsv

# protein files that have the functional annotations
outprefix=candidate
func_prot=${outprefix}_anno.faa
func_prot_Pdom_only=${outprefix}_Pdoms.anno.faa

evalue_cutoff="1e-6"
threads=32

log_file=functional_anno_analysis.log

unset in_run_anno

#### basic functions ########################################################################################################

function msg { # write a msg to stderr
   >&2 echo -e $@
}

function log {  # output msg and if in the run_anno process write to the log file
   msg $@
   [ -z "$in_run_anno" ] && return

   echo -e $(date +"%a %d-%b-%Y %r"): "$@" >> $log_file
}

function blank_log_line {
   msg ""
   [ -z "$in_run_anno" ] && return

   echo "" >> $log_file
}

function pprt_tm {  # input is in seconds
  awk -v in_secs=$1 '
    function pprt_tm(S    ,s,m,h) {
        s=(S%60);  tm = (s>9) ? s"s" : "0"s"s"
        if(S >= 60) {
            m=int(S/60)
            if(m >= 60) {
               h=int(m/60)
               if(h>=24) { d=int(h/24); h=d"d"(h % 24) }
               mnt=m % 60; if(mnt<10) mnt = "0"mnt
               m = h"h"mnt
            }
            tm = m"m" tm
        }
        return tm
    }
    BEGIN { print pprt_tm(in_secs) }
  '
}
function pprt_span { # input is start_time in seconds and end_time in seconds
   start_time=$1
   end_time=$2
   [ -z $end_time ] && return 1

   seconds=$(awk -v start_time=$start_time -v end_time=$end_time 'BEGIN{ print end_time-start_time }')
   pprt_tm $seconds
}
function systime {
   awk 'BEGIN{print systime()}'
}

function start_run {
   unset run_time
   start_time=$(systime)
}
function end_run {
   [ -z $start_time ] && unset run_time && return

   end_time=$(systime)
   run_time=$( awk -v start_time=$start_time -v end_time=$end_time 'BEGIN{ print end_time-start_time }' )

   unset start_time
}

: ' # InterProScan-5.57-90.0

Available analyses:
                      TIGRFAM (15.0) : TIGRFAMs are protein families based on hidden Markov models (HMMs).
                         SFLD (4) : SFLD is a database of protein families based on hidden Markov models (HMMs).
                  SUPERFAMILY (1.75) : SUPERFAMILY is a database of structural and functional annotations for all proteins and genomes.
                      PANTHER (15.0) : The PANTHER (Protein ANalysis THrough Evolutionary Relationships) Classification System is a unique resource that classifies genes by their functions, using published scientific experimental evidence and evolutionary relationships to predict function even in the absence of direct experimental evidence.
                       Gene3D (4.3.0) : Structural assignment for whole genes and genomes using the CATH domain structure database.
                        Hamap (2021_04) : High-quality Automated and Manual Annotation of Microbial Proteomes.
                        Coils (2.2.1) : Prediction of coiled coil regions in proteins.
              ProSiteProfiles (2022_01) : PROSITE consists of documentation entries describing protein domains, families and functional sites as well as associated patterns and profiles to identify them.
                        SMART (7.1) : SMART allows the identification and analysis of domain architectures based on hidden Markov models (HMMs).
                          CDD (3.18) : CDD predicts protein domains and families based on a collection of well-annotated multiple sequence alignment models.
                       PRINTS (42.0) : A compendium of protein fingerprints - a fingerprint is a group of conserved motifs used to characterise a protein family.
                        PIRSR (2021_05) : PIRSR is a database of protein families based on hidden Markov models (HMMs) and Site Rules.
              ProSitePatterns (2022_01) : PROSITE consists of documentation entries describing protein domains, families and functional sites as well as associated patterns and profiles to identify them.
                      AntiFam (7.0) : AntiFam is a resource of profile-HMMs designed to identify spurious protein predictions.
                         Pfam (35.0) : A large collection of protein families, each represented by multiple sequence alignments and hidden Markov models (HMMs).
                   MobiDBLite (2.0) : Prediction of intrinsically disordered regions in proteins.
                        PIRSF (3.10) : The PIRSF concept is used as a guiding principle to provide comprehensive and non-overlapping clustering of UniProtKB sequences into a hierarchical order to reflect their evolutionary relationships.

Deactivated analyses:
                      Phobius (1.01) : Analysis Phobius is deactivated, because the resources expected at the following paths do not exist: bin/phobius/1.01/phobius.pl
                  SignalP_EUK (4.1) : Analysis SignalP_EUK is deactivated, because the resources expected at the following paths do not exist: bin/signalp/4.1/signalp
        SignalP_GRAM_POSITIVE (4.1) : Analysis SignalP_GRAM_POSITIVE is deactivated, because the resources expected at the following paths do not exist: bin/signalp/4.1/signalp
        SignalP_GRAM_NEGATIVE (4.1) : Analysis SignalP_GRAM_NEGATIVE is deactivated, because the resources expected at the following paths do not exist: bin/signalp/4.1/signalp
                        TMHMM (2.0c) : Analysis TMHMM is deactivated, because the resources expected at the following paths do not exist: bin/tmhmm/2.0c/decodeanhmm, data/tmhmm/2.0c/TMHMM2.0c.model

Some applicable options for us
 -appl,--applications <ANALYSES>                           Optional, comma separated list of analyses.  If this option
                                                           is not set, ALL analyses will be run.
 -exclappl,--excl-applications <EXC-ANALYSES>              Optional, comma separated list of analyses you want to exclude.
 -dra,--disable-residue-annot                              Optional, excludes sites from the XML, JSON output
'

function run_interproscan {
   [ ! -d "ipr" ] && mkdir ipr && log "Created ipr directory for interproscan results"
   date > ${prefix}_start

   exclude_these="MobiDBLite,PIRSR,SFLD,Coils,Gene3D,Hamap"

   log "running interproscan"
   interproscan.sh -i $protein_file -b $prefix --excl-applications $exclude_these -cpu $threads --disable-residue-annot -goterms -f TSV -f XML -f GFF3 -verbose |& tee ${prefix}.log
   retVal=$?
   [ $retVal -ne 0 ] && return 1
   [ ! -s $ipr_tsv ] && return 1 # return if the tsv file does not exist or exists and is empty

   # looks like interproscan completed successfully
   date > $ipr_end_file

   # make a file that has DbxRef info gleaned from the Interproscan tsv file
   log "pulling the Dbxref info from $ipr_tsv and storing in $ipr_anno"
   iprparse.sh $ipr_tsv >$ipr_anno
}

function fasta_subset_by_Pfam_IPR {
   # using the tsv format interproscan results, only keep those fasta records in the input file of arg2
   # where there has been a Pfam match or another match that identifies an IPR id

   interproscan_tsv=$1
   fasta_to_subset=$2

   fold_width=90

   bawk '
      NR!=FNR && $name in keep {
         print ">"$name " " $comment
         print $seq
         next
      }

      /Pfam/{keep[$1]++}/IPR[0-9]*/{keep[$1]++}

    ' <(prefix_gt $interproscan_tsv) $fasta_to_subset | fold -w $fold_width
}

function tsv_stats {

   tsv=$1
   fld=$2
   [ -z $fld ] && fld=11

   [ ! -z $run_time ] && time_msg=" in $(pprt_tm $run_time)"
   [   -z $run_time ] && time_msg="."

   # print header with count of queries found and time it took
   stats_header_msg=$(awk -v fld=$fld -v time_msg="$time_msg" '
      $1==lst{next}{lst=$1; query++}
      END{print query, "queries found" time_msg}
   ' $tsv)
   echo -e "${stats_header_msg}\n\neValue distribution:"

   # print out stats for the tsv best hit
   # count 0's separately and then sort from smallest e-values to largest and print those counts

   awk -v fld=$fld '
      $1==lst{next}{lst=$1}
      $fld ~ "^0"{print $fld}' $tsv \
    | sort | uniq -c

   awk -v fld=$fld '
      $1==lst{next}{lst=$1}
      $fld ~ "^0"{next}
      {sub("^.*e-", "e-",$fld)}{print $fld}' $tsv \
    | sort | uniq -c | sort -k2,2Vr
}

function make_stats_file {
   tsv=$1; db=$2
   unset stats_header_msg

   bn=$(basename $tsv .tsv)
   stats_file=${tsv_dir}/${bn}.stats

   # 22Nov2022 handle empty tsv file
   if [ ! -s $tsv ]; then
      rm -f $stats_file
      log "WARNING: no hits found in $tsv"
      return 1
   fi

   # 26Jan2024 remember db info if db name passed in
   >$stats_file
   [ ! -z $db ] && blastdbinfo $db >$stats_file && echo "" >>$stats_file

   tsv_stats $tsv >>$stats_file

   [ ! -z "$stats_header_msg" ] && log "$tsv ${stats_header_msg}\n"
}

function stats_file_exists {
   tsv=$1

   bn=$(basename $tsv .tsv)
   stats_file=${bn}.stats

   [ -s $stats_file ] || [ -s $tsv_dir/$stats_file ]
}

function tsv_complete {  # we presume the tsv file is complete if the stats file exists
   tsv=$1
   stats_file_exists $tsv
}

# tallies species that are in brackets at end of the fasta annotation, as in candidate_Pdoms.anno.faa
function species_tally {
   annotated_fasta=$1
   grep "\[[^\[]*\]$" $annotated_fasta -o |sort |uniq -c |sort -nr
}

# call with database and outfile arguments
function blastp_db {

   db=$1
   outfile=$2
   start_run

   time blastp -query $protein_file                      \
            -db $db                                      \
            -out $outfile                                \
            -outfmt "6 std qlen slen staxids stitle"     \
            -evalue $evalue_cutoff -max_hsps 1 -max_target_seqs 5  \
            -num_threads $threads -mt_mode 1 -seg yes

   retVal=$?
   end_run
   [ $retVal -eq 0 ] && make_stats_file $outfile $db
}

# call with database and outfile arguments
function diamond_blastp_more-sensitive {
   # this version uses the --more-sensitive flag
   # running says: The host system is detected to have 2160 GB of RAM. It is recommended to use this parameter for better performance: -c1
   # so we have added the -c1 flag
   # stitle needs to be the 16th tabbed delimited field, and can not get taxids into diamond db for orthodb formats, so use qstrand as placeholder

   db=$1
   outfile=$2
   start_run

   std="qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

   time diamond blastp -v --query $protein_file            \
            --db $db                                       \
            --out $outfile                                 \
            --outfmt 6 $std qlen slen qstrand stitle       \
            --evalue $evalue_cutoff --max-hsps 1 --max-target-seqs 5 \
            --more-sensitive -c 1 --threads $threads

   retVal=$?
   end_run
   [ $retVal -eq 0 ] && make_stats_file $outfile $db
}

function blast_nt {
   log "blastn nt database"
   qry=$cds_file
   outfile=$tsv_dir/$nt_tsv
   start_run

   time blastn -db nt -query $qry -out $outfile                              \
                 -outfmt "6 std staxid stitle qlen qcovhsp qcovus"           \
                 -evalue $evalue_cutoff -max_target_seqs 5 -subject_besthit  \
                 -num_threads $threads

   retVal=$?
   end_run

   if [ $retVal -eq 0 ]; then
      make_stats_file $outfile nt
      categorized=${tsv_dir}/categorized_${nt_tsv}
      categorize_tsv.sh $outfile >$categorized
      log "create $categorized"
   fi
}

function get_fna_recs_unmatched_by_blastps {
   bawk '
      NR==FNR {
         matched[$1]++
         next
      }
      ! ($name in matched) {  # write it out, it was not seen in any of our blast outputs
         print ">"$name " " $comment
         print $seq
      }
   ' <(prefix_gt *blastp.tsv) $cds_file > $unmatched_cds_file
}

function blast_nt_unmatched_recs {
   qry=$unmatched_cds_file
   outfile=$tsv_dir/$unmatched_nt_tsv
   start_run

   time blastn -db nt -query $qry -out $outfile                              \
                 -outfmt "6 std staxid stitle qlen qcovhsp qcovus"           \
                 -evalue $evalue_cutoff -max_target_seqs 5 -subject_besthit  \
                 -num_threads $threads

   retVal=$?
   end_run

   if [ $retVal -eq 0 ]; then
      make_stats_file $outfile
      categorized=${tsv_dir}/categorized_${unmatched_nt_tsv}
      categorize_tsv.sh $outfile >$categorized
      log "create $categorized"
   fi
}

function get_faa_recs_unmatched_by_tsvs {
   bawk '
      NR==FNR {
         matched[$1]++
         next
      }
      ! ($name in matched) {  # write it out, it was not seen in any of our blast outputs
         print ">"$name " " $comment
         print $seq
      }
   ' <(prefix_gt *.tsv) $protein_file > $unmatched_protein_file
}
#############################################################################################################################

#### specific functions calling the basic functions #########################################################################
function run_and_report_interproscan {
   if [ ! -s $ipr_tsv ]; then

      run_interproscan

      if [ ! -s $ipr_tsv ]; then
         log "interproscan did not finish successfully (no, or empty, file $ipr_tsv)"
      elif [ ! -s $ipr_anno ]; then
         log "$ipr_tsv file exists but $ipr_anno does not: pull the Dbxref info from $ipr_tsv and store in $ipr_anno"
         iprparse.sh $ipr_tsv >$ipr_anno
      fi
   else
      if [ ! -s $ipr_anno ]; then
         log "$ipr_tsv file exists but $ipr_anno does not: pull the Dbxref info from $ipr_tsv and store in $ipr_anno"
         iprparse.sh $ipr_tsv >$ipr_anno
      else
         log "found \"$ipr_anno\", interproscan run complete."
      fi
   fi
}

function blastp_orthodb_vertebrates {
   log "blastp orthodb_vertebrates"
#   blastp_db /ccg/bin/orthodb_downloads/vertebrate_aa $tsv_dir/orthodb_verts.blastp.tsv
   blastp_db /ccg/bin/orthodb_downloads/Vertebrata_OrthoDB_11 $tsv_dir/orthodb_verts.blastp.tsv
}

function blastp_orthodb_arthropoda {
   log "blastp orthodb_arthropoda"
   blastp_db /ccg/bin/orthodb_downloads/arthropoda_orthodb $tsv_dir/orthodb_arthropoda.blastp.tsv
}
function diamond_blastp_orthodb_vertebrates {
   log "diamond blastp orthodb_vertebrates"
#   diamond_blastp_more-sensitive /ccg/bin/orthodb_downloads/vertebrate_aa $tsv_dir/orthodb_verts.blastp.tsv
   diamond_blastp_more-sensitive /ccg/bin/orthodb_downloads/Vertebrata_OrthoDB_11 $tsv_dir/orthodb_verts.blastp.tsv
}

function diamond_blastp_orthodb_arthropoda {
   log "diamond blastp orthodb_arthropoda"
   diamond_blastp_more-sensitive /ccg/bin/orthodb_downloads/arthropoda_orthodb $tsv_dir/orthodb_arthropoda.blastp.tsv
}

function blastp_swissprot {
   log "blastp SwissProt"
   blastp_db swissprot $tsv_dir/uniprot_sprot.blastp.tsv  # was uniprot_sprot
}

function blastp_trembl {
   log "blastp TrEMBL"
   blastp_db uniprot_trembl $tsv_dir/uniprot_trembl.blastp.tsv
}

function diamond_blastp_trembl {
   log "diamond blastp TrEMBL"
   diamond_blastp_more-sensitive /ccg/db_sets/uniprot/trembl $tsv_dir/uniprot_trembl.dmnd.more-sensitive_blastp.tsv
}

function diamond_blastp_nr {
   log "diamond blastp nr"
   diamond_blastp_more-sensitive /ccg/blastdbs/nr $tsv_dir/nr.dmnd.more-sensitive_blastp.tsv
}

function eggNOG_analysis {
   mkdir -p eggnog

   min_eval=0.000001

   cmd="eggnog.sh -m diamond --sensmode more-sensitive -i $protein_file -o eggnog/cogs --dbmem --seed_ortholog_evalue $min_eval --cpu $threads"
   log "eggNOG analysis: $cmd"

   $cmd
}

#############################################################################################################################

function set_orthodb {
   orthodb_tsv="_"
   if [ ! -z $vertebrate ]; then
      orthodb_tsv=$tsv_dir/orthodb_verts.blastp.tsv
      orthodb_blast_func=blastp_orthodb_vertebrates
   elif [ ! -z $arthropoda ]; then
      orthodb_tsv=$tsv_dir/orthodb_arthropoda.blastp.tsv
      orthodb_blast_func=blastp_orthodb_arthropoda
   fi
   if [ ! $orthodb_tsv == "_" ] && [ ! -z $diamond_orthodb ]; then
      orthodb_blast_func="diamond_$orthodb_blast_func"
   fi
}

function extract_recs_unmatched_by_tsvs_then_run_blastn_with_nt {
   log "extracting fna records unmatched by blastp runs"
   get_fna_recs_unmatched_by_blastps
   msg $unmatched_cds_file: $(numrecs $unmatched_cds_file) records

   log "blastn nt on the unmatched records"
   blast_nt_unmatched_recs
}
: '
function extract_recs_unmatched_by_tsvs_then_run_blastp_with_nr {
   log "extracting faa records unmatched by previous runs"
   get_fna_recs_unmatched_by_blastps
   msg $unmatched_cds_file: $(numrecs $unmatched_cds_file) records

   log "blastp nr on the unmatched records"
   blast_nt_unmatched_recs
}
'

function subset_Pfam_IPR_recs_for_search {
   new_prefix=PF_IPR_

   out_faa=${new_prefix}$protein_file
   fasta_subset_by_Pfam_IPR $ipr_tsv $protein_file >$out_faa
   protein_file=$out_faa

   out_cds=${new_prefix}$cds_file
   fasta_subset_by_Pfam_IPR $ipr_tsv $cds_file >$out_cds
   cds_file=$out_cds

   faa_recs=$(numrecs $protein_file)
   log "rest of analysis will only use those sequences with a Pfam domain or an Interproscan IPR id\n$orig_recs sequences originally\n$faa_recs with Pfam or IPR"
}

function write_ipr_distributions {
   # define the result gff as being the gff file most recently written to
   out_gff=$(ls -t *gff | head -1)
   [ -z $out_gff ] && log "Could not find output gff name to use for ipr distribution survey" && return 14

   dist_file=ipr_distributions.txt

   ipr_dist.sh $out_gff > $dist_file
   [ -s $dist_file ] && log "$dist_file created from mRNA lines in $out_gff"
}

function tally_annotated_faa {
   unset tally_faa
   tally_file=species_tally.txt

   [ -s $func_prot ] && tally_faa=$func_prot
   [ -s $func_prot_Pdom_only ] && tally_faa=$func_prot_Pdom_only

   [ ! -z $tally_faa ] && species_tally $tally_faa > $tally_file

   if [ -s $tally_file ]; then
      log "$tally_file created from $tally_faa annotations"
   else
      log "$tally_file was not created"
   fi
}

function basic_stats {
   # define the result gff as being the gff file most recently written to
   out_gff=$(ls -t *gff | head -1)
   [ -z $out_gff ] && log "Could not find output gff name to use for basic stats" && return 15

   pct_msg="; to add percentages you can run with assembly length as shown here: basic_gff_stats.sh $out_gff <assembly file>"
   [ -s assembly ] && pct_msg=""

   basic_stat_file=basic_gff_stats.txt
   basic_gff_stats.sh $out_gff assembly >$basic_stat_file
   [ -s $basic_stat_file ] && log "$basic_stat_file has basic gene and exon stats${pct_msg}"
}

function create_info_files {
   write_ipr_distributions # maybe useful count of the various tools
   tally_annotated_faa
   : basic_stats  # done for each of the output gffs in make_gene_model_tsv_and_anno_gff.sh now, do not need it here
}

function make_gene_models {
   log "running make_gene_model_tsv_and_anno_gff.sh with nt SwissProt TrEMBL OrthoDB InterProscan tsvs"

   gm_dir=$tsv_dir
   set_orthodb

   sprot_tsv=$gm_dir/uniprot_sprot.blastp.tsv
   trembl_tsv=$gm_dir/uniprot_trembl.dmnd.more-sensitive_blastp.tsv
   nt_tsv=$gm_dir/$nt_tsv

   mgm_log=$log_file
   [ -z "$in_run_anno" ] && mgm_log="/dev/null"
   log "   make_gene_model_tsv_and_anno_gff.sh $nt_tsv $sprot_tsv $trembl_tsv $orthodb_tsv $ipr_anno"

   # NT=$1; Sprot=$2; TrEMBL=$3; OrthoDB=$4; IPR=$5
   make_gene_model_tsv_and_anno_gff.sh $nt_tsv $sprot_tsv $trembl_tsv $orthodb_tsv $ipr_anno | tee -a $mgm_log
}
#############################################################################################################################

function run_blast_create_tsv {
   func_to_call=$1
   file_to_write=$2

   if ! tsv_complete $file_to_write; then
      $func_to_call
   else
      log "$file_to_write already complete"
   fi
}

### main routine to call ####################################################################################################

function run_anno {

   [ -z $threads ] && threads=32

   start_anno_run=$(systime)
   start_msg="Start functional annotation analysis"
   [ -s $log_file ] && start_msg="Restarting annotation run"

   in_run_anno="$start_msg"

   log $start_msg

   # handle interproscan run
   run_and_report_interproscan

   orig_recs=$(numrecs $protein_file)

   # if pf_ipr_only set, do blast analyses only on the PF/IPR subset (usually best to do on all and exclude non matches in the make_gene_models step)
   if [ ! -z $pf_ipr_only ] && [ -f $ipr_tsv ]; then
      subset_Pfam_IPR_recs_for_search
   else
      faa_recs=$orig_recs
   fi

   # functions use $protein_file and $cds_file in the searches

   [ ! -d $tsv_dir ] && mkdir $tsv_dir && log "Created $tsv_dir directory for blastn, blastp and diamond results"

   # 1 nt
   run_blast_create_tsv blast_nt $nt_tsv

   # 2 eggnog   14Oct2023 add eggNOG analysis
   [ ! -s eggnog/cogs.emapper.annotations ] && eggNOG_analysis

   # 3 swissprot
   run_blast_create_tsv blastp_swissprot uniprot_sprot.blastp.tsv

   #4 TrEMBL
   run_blast_create_tsv diamond_blastp_trembl uniprot_trembl.dmnd.more-sensitive_blastp.tsv

   # 5 orthodb
   set_orthodb  # set vars for orthodb processing
   if [ $orthodb_tsv == "_" ]; then
      log "no orthodb run requested"
   elif tsv_complete $orthodb_tsv; then
      log "$(basename $orthodb_tsv) already complete"
   else
      eval $orthodb_blast_func # blastp_orthodb_vertebrates or blastp_orthodb_arthropoda
   fi

   # 6 nr  14Oct2023 add back NR search. it takes 10 hours or more but adds info not found in other searches
   run_blast_create_tsv diamond_blastp_nr nr.dmnd.more-sensitive_blastp.tsv

   # analyze blasts and make gene-models, gffs, faa and fna files
   make_gene_models

   create_info_files

   blank_log_line
   run_time=$(pprt_span $start_anno_run $(systime))

   log "End of functional annotation run. $run_time\n"
}
#############################################################################################################################
