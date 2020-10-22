#!/bin/bash

################################################################################
# Info
################################################################################

# SARSeq analysis script to count known amplicons per plate and well from an
# Illumina paired-end (PE) run with 2 additional index reads (4 fastq.gz files) 
#
# The script reports raw read counts and does NOT predict or diagnose infection
# or disease, which requires an expert interpreting the read counts for viral 
# and control amplicons (similar to RT-qPCR Ct values that also require expert
# interpretation).
#
# In our experience, negative samples have 0 or a few counts only, while positive
# samples have hundreds or thousands of counts depending on virus titer and sequencing
# depth.
#
# See Yelagandula et al., bioRxiv 2020 for details on SARSeq and all results.

# Based on SARSeq's two-dimensional unique dual indexing and primer design,
# Plate ID is redundantly encoded by the i5 and i7 Illumina indices (2nd PCR).
# Well ID is reduntantly encoded by indices in the forward and reverse primers
# of the 1st PCR, whereby the forward well index is after a staggered offset of 1-4 nts
# this offset can either be random or fixed for all amplicons of a particular well
#
# forward primer: 5'adapter-offset-wellID-forward
# reverse primer: 3'reverse-wellID-adapter
#
# forward read: 5'offset-wellID-amplicon (wellID starts at positions 2 to 5)
# reverse read: 5'wellID-amplicon (no random offset)
#
# Developed for ease of use on any Unix/Linux system, with minimal dependencies.
# It therefore conciously does not make use of common tools in genomics and NGS
# data analysis and allows at most 1 mismatch per sequence (i.e. each of the 4
# indices and the amplicon). In our tests, 0 mismatches for the amplicon gave the
# highest specificity.
#
#
# Written in summer/fall 2020 by
#
# Alex Stark
# Research institute of Molecular Pathology (IMP)
# Vienna, Austria
# 
# Alex thanks Vanja Haberle, Bernardo Almeida and Lisa-Marie Pleyer (all IMP)
# for testing the script
#
# Released under the MIT license


################################################################################
# Requirements
################################################################################

# System:
#
# Linux/Unix OS
# zcat/gunzip
#
# 6 cores recommended

# Files:
#
# i5i7_NextSeq.txt or i5i7_MiSeq.txt (plate IDs)
# wellIDs.txt or wellIDs_fixedOffsets.txt
# amplicons.fa or other_viruses.fa or all_amplicons.fa


################################################################################
# Set default values
################################################################################

mm=0         # number of mismatches in amplicon (0/1)
MM=0         # number of mismatches in plate and well indices (0/1)

LEN=0       # length of amplicon to be considered (shortest amplicon minus 3)
            # depending on provided amplicon fasta file: 69 or 55
            # default 0 means that it will be determined from the fasta file


SCRIPTDIR=$(dirname $(readlink -f $0)) # path to this script, where also the indices and fasta files are found

IDX=""                                 # text file with i5/i7-to-plate mappings (i5i7_NextSeq.txt or raw/i5i7_MiSeq.txt)
WLL=$SCRIPTDIR/wellIDs.txt             # text file with well IDs
FIXEDOFF="0"                           # offset is fixed (0/1); fixed offset requires length of offset (integer) in col 2 of WLL

plateMin=101 # min plate index to report
plateMax=496 # max plate index to report

FA=$SCRIPTDIR/amplicons.fa # fasta file of amplicons to use

DIR=$(pwd) # make the default the current working directory

NOCOPY=0 # hidden parameter: don't copy any fastq files, just assume that the fastq files are already there


################################################################################
# Help
################################################################################

print_help () {
    echo >&2 "

$(basename $0) - SARSeq analysis script

USAGE: $(basename $0) -1 <forward-reads.fastq.gz> -2 <reverse-reads.fastq.gz> -3 i5-reads.fastq.gz -4 i7-reads.fastq.gz -I i5/i7-indices -t runID [OPTIONS] 

 -1     fastq.gz file forward reads (typically named Undetermined..._R1_001.fastq.gz) [required]
 -2     fastq.gz file reverse reads (typically named Undetermined..._R2_001.fastq.gz) [required]
 -3     i5 indices as separate fastq.gz file (typically named Undetermined..._I1_001.fastq.gz) [required]
 -4     i7 indices as separate fastq.gz file (typically named Undetermined..._I2_001.fastq.gz) [required]

 -I     text file with i5/i7 index-to-plate mapping [required]
 -t     run ID (alphanumerical string; 3-digit ID 001 to 999 recommended) [required]

 -m     mismatches in amplicons (0 or 1) [default: $mm]
 -M     mismatches in plate and well indices (0 or 1) [default: $MM]
 -l     length of amplicon considered (shortest amplicon [minus 3 for random offset]; 0: autom. determined) [default: $LEN]

 -i     min plate index to report [default: $plateMin]
 -a     max plate index to report [default: $plateMax]

 -A     fasta file with amplicons [default: $FA]
 -W     text file with well indices [default: $WLL] 
 -F     offset is fixed for all amplicons per well (binary: 0/1) [default: $FIXEDOFF]
        -F 1 allows the reporting of the number of unmatched and total reads per well/sample
        and in our experience slightly increases the specificty of read-to-well mapping 

 -d     output directory [default: $DIR]


NOTES:

This script counts SARSeq reads from a paired-end NGS run for each sample/well and amplicon, reporting raw read counts.

It does NOT predict or diagnose infection or disease - this requires an expert interpreting the read counts for viral 
and control amplicons, similar to RT-qPCR Ct values that also require expert interpretation.

Choose i5/i7 index file with -I as $SCRIPTDIR/i5i7_NextSeq.txt or $SCRIPTDIR/i5i7_MiSeq.txt depending on sequencing run type

Recommended use of 6 compute cores

Output files will be generated in $DIR/[runID]/ - this directory must not exist

Two types of output files will be generated:
- plates_*.txt -> readcounts in 96-well format; one file per amplicon
- table_*.txt  -> readcounts in table format; one single file, one sample/well per row, one amplicon per column

"

}


################################################################################
# Parse input
################################################################################

if [ $# -eq 0 ]; then
    print_help;  
    exit 1
fi

while getopts "1:2:3:4:t:m:M:l:i:a:A:I:W:F:d:Ch" o
do
  case "$o" in
    1) FQF="$OPTARG";;
    2) FQR="$OPTARG";;
    3) FQi5="$OPTARG";;
    4) FQi7="$OPTARG";;
    I) IDX="$OPTARG";;
    t) RN="$OPTARG";;
    m) mm="$OPTARG";;
    M) MM="$OPTARG";;
    l) LEN="$OPTARG";;
    i) plateMin="$OPTARG";;
    a) plateMax="$OPTARG";;
    A) FA="$OPTARG";;
    W) WLL="$OPTARG";;
    F) FIXEDOFF="$OPTARG";;
    d) DIR="$OPTARG";;
    C) NOCOPY=1;;
    h) print_help; exit;;
   \?) print_help; exit 1;;
  esac
done

## check if index file is given and exists
if [ -z "$IDX" ]; then
    echo >&2 "ERROR: -I [i5/i7 index file] is required!  
Choose between $SCRIPTDIR/i5i7_NextSeq.txt or $SCRIPTDIR/i5i7_MiSeq.txt depending on sequencing run type"
    exit 1
elif [ ! -f "$IDX" ]; then
    echo >&2 "ERROR: i5/i7 index file $IDX does not exist!
Choose between $SCRIPTDIR/i5i7_NextSeq.txt or $SCRIPTDIR/i5i7_MiSeq.txt depending on sequencing run type"
    exit 1
fi

## check if run ID is given and has only alphanumerical characters
if [ -z "$RN" ]; then
  echo >&2 "ERROR: -t [run ID] is required!"; exit 1
fi
NC=$(awk -v RN=$RN 'BEGIN{if(RN~/^[0-9a-zA-Z]+$/){print 1}else{print 0}}')
if [ "$NC" -ne "1" ]; then
    echo >&2 "ERROR: -t [run ID] must only contain alphanumeric characters!"; exit 1
fi

## check WLL file exists and is in the correct format given the offset choice
if [ ! -f "$WLL" ]; then
    echo >&2 "ERROR: well index file $WLL does not exist!"; exit 1
fi
NC=$(awk -F'\t' '{if(NR==1){N=NF}else{if(N!=NF){N=-1}}} END{print N}' $WLL)
if [ "$FIXEDOFF" -eq "0" ]; then
    # this is the version without fixed offset
    if [ "$NC" -ne "3" ]; then
        echo >&2 "ERROR: wrong format for well index file $WLL
When running without fixed offset, 3 cols are required: well-ID, fw-index, rv-index"
        exit 1
    fi
else 
    # this is the version with fixed offset
    if [ "$NC" -ne "4" ]; then
        echo >&2 "ERROR: wrong format for well index file $WLL
When running with fixed offset, 4 cols are required: well-ID, offset, fw-index, rv-index"
        exit 1
    fi
fi

## check if mm and MM parameters are sensible
if [ $(awk -v mm="$mm" -v MM="$MM" 'BEGIN{N[0]=1;N[1]=1;N["N"]=1;if((mm in N) && (MM in N)){print 1}else{print 0}}') -eq 0 ]; then 
  echo >&2 "ERROR: mm and MM must be 0 or 1 (but are $mm and $MM)"; exit 1
fi

## check if LEN is given, otherwise determine from fasta file
if [ "$LEN" -le "0" ]; then
    LEN=$(awk -vFA=$FA -vFO=$FIXEDOFF 'BEGIN{while((getline<FA)>0){if(substr($1,1,1)==">"){id=substr($1,2);LEN[id]=0}else{LEN[id]=LEN[id]+length($1)}}
                                             MINLEN=1000000;for(id in LEN){if(LEN[id]<MINLEN){MINLEN=LEN[id]}} if(FO==1){print MINLEN}else{print MINLEN-3}}')
fi

## check hidden parameter and existance of fastq.gz files
if [ "$NOCOPY" -eq "0" ]; then
  if [ -z "$FQF" -o -z "$FQR" -o -z "$FQi5" -o -z "$FQi7" ]; then
    echo >&2 "ERROR: fastq.gz input files -1 -2 -3 -4 are required!"; exit 1
  fi
  if [ ! -f "$FQF" -o ! -f "$FQR" -o ! -f "$FQi5" -o ! -f "$FQi7" ]; then 
    echo >&2 "ERROR: input files $FQF, $FQR, $FQi5 or $FQi7 not found"; exit 1 
  fi 
  if [ -e "$DIR/$RN" ]; then
    echo >&2 "ERROR: $DIR/$RN must not exist"; exit 1
  fi
fi


################################################################################
# Start analysis
################################################################################

## create analysis folder
mkdir -p "$DIR/$RN/raw"

echo "starting analysis with parameters
     RUN: $RN
     mm: $mm
     MM: $MM 
     LEN: $LEN
     FA: $FA
     IDX: $IDX
     WLL: $WLL

     FixedOffset: $FIXEDOFF

     ...

     ... creating analysis folder $DIR/$RN" > $DIR/$RN/log.${mm}.${MM}

## cp fastq files to analysis folder
if [ "$NOCOPY" -eq "0" ]; then
  echo "... copying fastq files to analysis folder $DIR/$RN/raw 
            copying $FQF
            copying $FQR 
            copying $FQi5
            copying $FQi7" >> $DIR/$RN/log.${mm}.${MM}
  PIDS=""
  cp $FQF $DIR/$RN/raw/Undetermined_S0_R1_001.fastq.gz &
  PIDS="$PIDS $! "
  cp $FQR $DIR/$RN/raw/Undetermined_S0_R2_001.fastq.gz &
  PIDS="$PIDS $! "
  if [ ! -z "$FQi5" -a -f "$FQi5" ]; then
    cp $FQi5 $DIR/$RN/raw/Undetermined_S0_I1_001.fastq.gz &
    PIDS="$PIDS $! "
  fi
  if [ ! -z "$FQi7" -a -f "$FQi7" ]; then
    cp $FQi7 $DIR/$RN/raw/Undetermined_S0_I2_001.fastq.gz &
    PIDS="$PIDS $! "
  fi
  wait $PIDS
else
  echo "... fastq files expected to be in analysis folder $DIR/$RN/raw" >> $DIR/$RN/log.${mm}.${MM}
fi

## define shorthand
FQ=$DIR/$RN/raw/Undetermined_S0

## go through fastq.gz files
echo "... mapping reads" >> $DIR/$RN/log.${mm}.${MM}

########## fork depending on FIXEDOFF

if [ "$FIXEDOFF" -eq "0" ]; then

##### this is the version without fixed offset
# Analsis strategy: use i5/i7 to get plate; get wellID once from reverse read
#                   determine amplicon and its position within forward read
#                   from this, determine random offset, and well barcode
#

## go through all 4 fastq.gz files in parallel    
paste <(zcat ${FQ}_R1_001.fastq.gz) <(zcat ${FQ}_R2_001.fastq.gz) <(zcat ${FQ}_I1_001.fastq.gz) <(zcat ${FQ}_I2_001.fastq.gz) | \
    awk -vOFS='\t' -vLEN=$LEN -vmm=$mm -vMM=$MM -vRN="$RN" -vFA=$FA -vIDX=$IDX -vWLL=$WLL 'BEGIN{
        NT["N"]=1;NT["A"]=1;NT["C"]=1;NT["G"]=1;NT["T"]=1
        while((getline<FA)>0){if(substr($1,1,1)==">"){ID=substr($1,2)}else{S[ID]=(S[ID] toupper($1))}}
        for(s in S){ for(o=1;o<=4;o++){ AMP[substr(S[s],o,LEN)]=s; OFF[substr(S[s],o,LEN)]=o } }
        if(mm>0){
            for(A in AMP){ for(i=1;i<=length(A);i++){ for(nt in NT){ AMP[(substr(A,1,i-1) nt substr(A,i+1))]=AMP[A]; OFF[(substr(A,1,i-1) nt substr(A,i+1))]=OFF[A] }}}
        }
        while((getline<IDX)>0){
            W5[$2]=$1;if(MM>0){for(i=1;i<=length($2);i++){for(nt in NT){W5[(substr($2,1,i-1) nt substr($2,i+1))]=$1}}}
            W7[$3]=$1;if(MM>0){for(i=1;i<=length($3);i++){for(nt in NT){W7[(substr($3,1,i-1) nt substr($3,i+1))]=$1}}} }
        while((getline<WLL)>0){
            Wf[$2]=$1;if(MM>0){for(i=1;i<=length($2);i++){for(nt in NT){Wf[(substr($2,1,i-1) nt substr($2,i+1))]=$1}}}
            Wr[$3]=$1;if(MM>0){for(i=1;i<=length($3);i++){for(nt in NT){Wr[(substr($3,1,i-1) nt substr($3,i+1))]=$1}}} }
    }
    ((NR%4==2) && (substr($1,13,LEN) in AMP)){
        w5=substr($4,1,8); w7=substr($3,1,8); wr=substr($2,1,8)
        if((w5 in W5) && (w7 in W7) && (W5[w5]==W7[w7]) && (wr in Wr)){
            amp=substr($1,13,LEN); wf=substr($1,6-OFF[amp],8)
            if((wf in Wf) && (Wf[wf]==Wr[wr])){
                N[(W5[w5]"-"Wf[wf]"\t"AMP[amp])]++
            }
        }
    }
    END{for(n in N){print (RN "-" n), N[n]}}' > $DIR/$RN/awk_mappings_${mm}.${MM}.txt

# check if mapping had results
if [ ! -s $DIR/$RN/awk_mappings_${mm}.${MM}.txt ]; then
    echo >&2 "Oops, something went wrong: mapping file $DIR/$RN/awk_mappings_${mm}.${MM}.txt is empty
Check if i5/i7 index file and other input files are correct"
    exit 1
fi

# results in plate format
awk '(/^>/){print substr($1,2)}' $FA | while read A; do
    cat $DIR/$RN/awk_mappings_${mm}.${MM}.txt | \
        awk -vA=$A -vPMin=$plateMin -vPMax=$plateMax -vRN=$RN '($2==A){N[$1]+=$3} 
                    END{for(p=PMin;p<=PMax;p++){print (">" p); L=""; for(c=1;c<=12;c++){if(c<10){C=("0"c)}else{C=c} L=(L "\t" C)} print L; 
                            for(r=1;r<=8;r++){R=substr("ABCDEFGH",r,1);L=R; 
                                for(c=1;c<=12;c++){if(c<10){C=("0"c)}else{C=c} L=(L"\t"N[RN"-"p"-"R C]+0)} print L}}}' > $DIR/$RN/plates_${A}_${mm}.${MM}.txt
done

# results in table format
awk -vAM=$DIR/$RN/awk_mappings_${mm}.${MM}.txt -vPMin=$plateMin -vPMax=$plateMax -vRN=$RN -F'\t' -vOFS='\t' \
    'BEGIN{while((getline<AM)>0){N[$1":"$2]=$3}; n=0} 
     (/^>/){n++;A[n]=substr($1,2)} 
     END{L="#RUN-PLATE-WELL";for(i=1;i<=n;i++){L=(L"\t"A[i])} print L
         for(p=PMin;p<=PMax;p++){
             for(r=1;r<=8;r++){R=substr("ABCDEFGH",r,1);
                 for(c=1;c<=12;c++){if(c<10){C=("0"c)}else{C=c} id=(RN"-"p"-"R C);L=id;for(i=1;i<=n;i++){L=(L"\t"N[id":"A[i]]+0)} print L}
             }
         }
     }' $FA > $DIR/$RN/table_${mm}.${MM}.txt


else
   
##### this is the version WITH fixed offset (from col 2 of WLL index)
# Analsis strategy: use i5/i7 to get plate; get wellID 1 from reverse read
#                   then determine fixed offset for that wellID
#                   then extract wellID 2 and amplicon via the known fixed offset
#                   determine amplicon identity (or "-")
# Note: this strategy can determine unmapped and thus "all reads" for each well

## go through all 4 fastq.gz files in parallel
paste <(zcat ${FQ}_R1_001.fastq.gz) <(zcat ${FQ}_R2_001.fastq.gz) <(zcat ${FQ}_I1_001.fastq.gz) <(zcat ${FQ}_I2_001.fastq.gz) | \
    awk -vOFS='\t' -vLEN=$LEN -vmm=$mm -vMM=$MM -vRN="$RN" -vFA=$FA -vIDX=$IDX -vWLL=$WLL 'BEGIN{
        NT["N"]=1;NT["A"]=1;NT["C"]=1;NT["G"]=1;NT["T"]=1
        while((getline<FA)>0){if(substr($1,1,1)==">"){ID=substr($1,2)}else{S[ID]=(S[ID] toupper($1))}}
        for(s in S){AMP[substr(S[s],1,LEN)]=s}
        if(mm>0){ 
            for(A in AMP){ for(i=1;i<=length(A);i++){ for(nt in NT){ AMP[(substr(A,1,i-1) nt substr(A,i+1))]=AMP[A] }}}
        }
        while((getline<IDX)>0){
            W5[$2]=$1;if(MM>0){for(i=1;i<=length($2);i++){for(nt in NT){W5[(substr($2,1,i-1) nt substr($2,i+1))]=$1}}}
            W7[$3]=$1;if(MM>0){for(i=1;i<=length($3);i++){for(nt in NT){W7[(substr($3,1,i-1) nt substr($3,i+1))]=$1}}} }
        while((getline<WLL)>0){
            O[$1]=$2+1
            Wf[$3]=$1;if(MM>0){for(i=1;i<=length($3);i++){for(nt in NT){Wf[(substr($3,1,i-1) nt substr($3,i+1))]=$1}}}
            Wr[$4]=$1;if(MM>0){for(i=1;i<=length($4);i++){for(nt in NT){Wr[(substr($4,1,i-1) nt substr($4,i+1))]=$1}}} }
    }
    (NR%4==2){
        w5=substr($4,1,8); w7=substr($3,1,8); wr=substr($2,1,8)
        if((w5 in W5) && (w7 in W7) && (W5[w5]==W7[w7]) && (wr in Wr)){
            wf=substr($1,O[Wr[wr]],8)
            if((wf in Wf) && (Wf[wf]==Wr[wr])){
                A=substr($1,O[Wr[wr]]+8,LEN);if(A in AMP){id=AMP[A]}else{id="-"}
                N[(W5[w5]"-"Wf[wf]"\t"id)]++
            }
        }
    }
    END{for(n in N){print (RN "-" n), N[n]}}' > $DIR/$RN/awk_mappings_${mm}.${MM}.txt


# check if mapping had results
if [ ! -s $DIR/$RN/awk_mappings_${mm}.${MM}.txt ]; then
    echo >&2 "Oops, something went wrong: mapping file $DIR/$RN/awk_mappings_${mm}.${MM}.txt is empty
Check if i5/i7 index file and other input files are correct"
    exit 1
fi

# results in plate format
awk 'BEGIN{print "AllReads\n-"} (/^>/){print substr($1,2)}' $FA | while read A; do
    cat $DIR/$RN/awk_mappings_${mm}.${MM}.txt | \
        awk -vA=$A -vPMin=$plateMin -vPMax=$plateMax -vRN=$RN '((A=="AllReads") || ($2==A)){N[$1]+=$3} 
                    END{for(p=PMin;p<=PMax;p++){print (">" p); L=""; for(c=1;c<=12;c++){if(c<10){C=("0"c)}else{C=c} L=(L "\t" C)} print L; 
                            for(r=1;r<=8;r++){R=substr("ABCDEFGH",r,1);L=R; 
                                for(c=1;c<=12;c++){if(c<10){C=("0"c)}else{C=c} L=(L"\t"N[RN"-"p"-"R C]+0)} print L}}}' > $DIR/$RN/${A}_plates_${mm}.${MM}.txt
done

# results in table format
awk -vAM=$DIR/$RN/awk_mappings_${mm}.${MM}.txt -vPMin=$plateMin -vPMax=$plateMax -vRN=$RN -F'\t' -vOFS='\t' \
    'BEGIN{while((getline<AM)>0){N[$1":"$2]=$3;N[$1":AllReads"]+=$3}; n=0} 
     (/^>/){n++;A[n]=substr($1,2)} 
     END{n++;A[n]="-";n++;A[n]="AllReads";L="#RUN-PLATE-WELL";for(i=1;i<=n;i++){L=(L"\t"A[i])} print L
         for(p=PMin;p<=PMax;p++){
             for(r=1;r<=8;r++){R=substr("ABCDEFGH",r,1);
                 for(c=1;c<=12;c++){if(c<10){C=("0"c)}else{C=c} id=(RN"-"p"-"R C);L=id;for(i=1;i<=n;i++){L=(L"\t"N[id":"A[i]]+0)} print L}
             }
         }
     }' $FA > $DIR/$RN/table_${mm}.${MM}.txt


fi ## end of FIXEDOFF fork

exit

