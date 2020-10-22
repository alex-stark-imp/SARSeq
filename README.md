SARSeq analysis script to count known amplicons per plate and well from an Illumina paired-end (PE) run with 2 additional index reads (4 fastq.gz files).

The script reports raw read counts and does NOT predict or diagnose infection or disease, which requires an expert interpreting the read counts for viral and control amplicons (similar to RT-qPCR Ct values that also require expert interpretation).

In our experience, negative samples have 0 or a few counts only, while positive samples have hundreds or thousands of counts depending on virus titer and sequencing depth.

See Yelagandula et al., bioRxiv 2020 for details on SARSeq and all results.

Based on SARSeq's two-dimensional unique dual indexing and primer design, Plate ID is redundantly encoded by the i5 and i7 Illumina indices (2nd PCR). Well ID is reduntantly encoded by indices in the forward and reverse primers of the 1st PCR, whereby the forward well index is after a staggered offset of 1-4 nts. This offset can either be random or fixed for all amplicons of a particular well

forward primer: 5'adapter-offset-wellID-forward
reverse primer: 3'reverse-wellID-adapter

forward read: 5'offset-wellID-amplicon (wellID starts at positions 2 to 5)
reverse read: 5'wellID-amplicon (no random offset)

Developed for ease of use on any Unix/Linux system, with minimal dependencies. It therefore conciously does not make use of common tools in genomics and NGS data analysis and allows at most 1 mismatch per sequence (i.e. each of the 4 indices and the amplicon). In our tests, 0 mismatches for the amplicon gave the highest specificity.

Written in summer/fall 2020 by

Alex Stark
Research institute of Molecular Pathology (IMP)
Vienna, Austria

Alex thanks Vanja Haberle, Bernardo Almeida and Lisa-Marie Pleyer (all IMP) for testing the script

Released under the MIT license
