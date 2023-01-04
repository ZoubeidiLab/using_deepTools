#!/bin/bash
#SBATCH --job-name="nrmChIP-NoCtrl"
#SBATCH -p express,big-mem,normal,long
#SBATCH --ntasks=1 --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=6:00:00
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --error=tmp/Job.%J.norm_ChIPseq.err
#SBATCH --output=tmp/Job.%J.norm_ChIPseq.out


###################################################################
###  Script written by Joshua M. Scurll, Vancouver Prostate Centre.
###################################################################


## Number of processors to use. Note that at least 8 processors are required.

nproc=32


## Name of Conda environment. deeptools and bedtools must BOTH be
## installed in this conda environment.

conda_env="deepTools_py3.8"



#########################################################
###  PATHS TO INPUT AND OUTPUT FILES AND DIRECTORIES  ###
#########################################################

## Path to a sorted ChIP-seq BAM file that does not have a corresponding
## input or IgG control.

bam1=""

## Path of directory for output bigWig files.
## DO NOT include the trailing slash (/) at the end of the path.

outdir="<...>/Normalized_bigWigs"



########################################
###  Options for RPGC normalization  ###
########################################

## BED file of regions to exclude from the computation of the scale factor for
## an average of 1x genomic coverage in the bam file. It is a good idea to
## exclude peak regions if peaks have already been identified. To not exclude
## any regions from RPGC computations, set RPGC_blacklist="None".

RPGC_blacklist="/groups/zoubeidigrp/downloads/genomic_blacklists/ENCODE/hg38_blacklist.ENCFF356LFX.sorted.merged.bed"

## Chromosomes to ignore during calculation of the RPGC normalization scale factor.
## If you skip any chromosomes here, then you will also need to adjust the effective
## genome size below accordingly (i.e. subtract the mappable size of the skipped chr).
## Any chromosomes specified here should be specified in the form of a comma-
## separated list of chromosome names with each chromosome quoted individually
## and then the entire list also quoted. Even if just one chromosome is given
## here, it should be nested within TWO sets of quotation marks.
##
## E.g. RPGC_chrsToIgnore="'chrX', 'chrY', 'chrM'" or SES_chrsToSkip="'chrM'"
##
## To keep all chromosomes, set this equal to the empty string "" or ''. Note
## that any chromosomes listed here should match the names of chromosomes in
## the bam files.

RPGC_chrsToIgnore="'chrM'"

## Effective genome size. This should be the total length of the genome that is
## mappable based on the read lengths and whether multi-mapping reads were removed
## (removing non-uniquely mapped regions of the genome reduces the effective genome
## size). If any blacklisted regions are excluded from the calculation of the RPGC
## scale factor, then the total length of all of these excluded regions should also
## be subtracted from the effective genome size, except when blacklisted regions
## overlap unmappable regions that have already been accounted for in the effective
## genome size.
##
## The code below automatically subtracts the total sum of blacklist region lengths
## from the value of eff_genome_size specified here by the user, so there is no need
## to do the subtraction manually. Note that this automatic calculation is NOT
## appropriate if the blacklisted regions include regions of low-mappability that
## have already been factored into the effective genome size. The size of ignored
## chromosomes will NOT be automatically subtracted.

eff_genome_size=$((2747877777-16569))  ## 16569 is the size of chrM.

if [ "${RPGC_blacklist}" != "None" ] && [ -f "${RPGC_blacklist}" ]; then
    blacklist_total_length=$( awk '{s+=($3-$2)}END{print s}' "${RPGC_blacklist}" )
    eff_genome_size=$(( eff_genome_size - blacklist_total_length ))
    echo "Total length of blacklisted regions has been subtracted from the user-specified"
    echo "    effective genome size. The new effective genome size is ${eff_genome_size}."
fi



#################################
###  Options for bamCoverage  ###
#################################

## Bin size:
bin_size=10

## Length for smoothing the ChIP track:
smooth_length=50

## Lengths for smoothing the ChIP track in order to generate the model for
## the background/biases. Specify this as a bash array of integers, e.g.
## input_smooth_lengths=(1010 10010).

input_smooth_lengths=(10010 100010)

## Relative weights for each value in input_smooth_lengths, given as a bash array.
## The background/bias model will be constructed as a weighted mean of different
## smoothed versions of the Input signal, and the relative weight of contribution
## of each of the smoothed Input tracks, in the order listed above, will be given
## by the corresponding value that you list here in the same order.
## E.g. if input_smooth_lengths=(110 1010 10010) and input_smooth_weights=(3 2 1),
## then the background will be ((3 x input smoothed over 110 bp) + (2 x input
## smoothed over 1010 bp) + (1 x input smoothed over 10010 bp)) / (3 + 2 + 1).

input_smooth_weights=(1 1)

## Extend reads? Set this to the length to which to extend single-end reads (or
## unmated reads in paired-end data). For paired end data, paired read mates will
## automatically be extended to the fragment length. To not extend reads, set this
## to "False".

read_ext=200

## Pseudocount for bigwigCompare ratio calculation. The pseudocount value(s) specified
## here will be added to the coverage value in each bin for bam1 and bam2 before
## computing the ratio in order to avoid divisions by 0. Note that the Input control
## bigWig file will be normalized to an average of 1x genome coverage (RPGC). Therefore,
## the pseudocount should be chosen to be somewhat smaller than 1. A value between 0.01
## and 0.5 (i.e. between 1% and 50% of the average Input coverage) could be appropriate.

pseudocount=0.1

## Additional options for bamCoverage as string:

opts="--ignoreDuplicates"



######################################
###  Suffix for output file names  ###
######################################

## The read extension length, bin size, smoothing length, and pseudocount will
## automatically be added to the file name. Specify here a file name suffix
## (exlcuding extension) for anything else descriptive that you wish to include
## in the filename. Leave this as an empty string "" or '' to just use the default
## output file name created without any additional suffix.

out_suffix=""



#######################################
###  Paths to UCSC utilities/files  ###
#######################################

## Specify path to directory containing the UCSC utility bedGraphToBigWig.
## DO NOT include trailing slash (/) at end of directory path.

ucsc_dir="/home/jscurll/Apps/UCSC_utilities"

## Specify path to a file containing the sizes of chromosome in the relevant
## reference genome, as obtained using the UCSC utility fetchChromSizes.

chrom_sizes_path="/groups/zoubeidigrp/downloads/UCSC/hg38.chrom.sizes"



#########################################################################
#########################################################################
#########################################################################



source ~/.bashrc
wait
conda activate ${conda_env}
wait


## Create output directory if it doesn't already exist.
if [ ! -d "${outdir}" ]; then
    mkdir -p "${outdir}"
fi

## Make a copy of this script in the output directory.
tstamp=$(date +%s)
cp -p $0 "${outdir}/BioinfoCmd.norm_ChIPseq_bam2bigwig.${tstamp}.txt"
chmod 755 "${outdir}/BioinfoCmd.norm_ChIPseq_bam2bigwig.${tstamp}.txt"


## Create suffix for output file names.

if [ "${read_ext}" != False ]; then
   out_suffix="${out_suffix}.e${read_ext}.bgRPGC"
else
   out_suffix="${out_suffix}.bgRPGC"
fi

if [ ${bin_size} -gt 1 ]; then
    out_suffix="${out_suffix}.bs${bin_size}"
fi

if [ ${smooth_length} -gt ${bin_size} ]; then
    out_suffix="${out_suffix}.smth${smooth_length}"
fi


## Put an extra set of quotation marks around blacklist file paths for use
## in Python commands.

if [ "${RPGC_blacklist}" != "None" ]; then
    RPGC_blacklist="'${RPGC_blacklist}'"
fi


## Compute the scale factor needed for RPGC normalization of bam1.

echo "Computing the RPGC scale factor for bam1..."

t0=$(date +%s)
RPGC_sf=$( python -c "
from deeptools.getScaleFactor import get_scale_factor
class args:
    bam='${bam1}'
    extendReads=${read_ext}
    scaleFactor=1.0
    filterRNAstrand=None
    normalizeUsing='RPGC'
    blackListFileName=${RPGC_blacklist}
    effectiveGenomeSize=${eff_genome_size}
    numberOfProcessors=${nproc}
    ignoreForNormalization=[${RPGC_chrsToIgnore}]
    minMappingQuality=False
    samFlagInclude=False
    samFlagExclude=False
    minFragmentLength=False
    maxFragmentLength=False
    verbose=False
RPGC_sf = get_scale_factor(args(), None)
print('\n')
print(RPGC_sf)
" | tail -n 1 | awk '{print $NF}' )
wait
t1=$(date +%s)

echo "Finished in $((t1-t0)) s."
echo "RPGC scale factor = ${RPGC_sf}"


## Construct the names of smoothed pseudo-Input tracks for use later.

bg_prefix="${outdir}/$(basename --suffix='.bam' "${bam1}")"
if [ "${read_ext}" != False ]; then
    bg_suffix=".e${read_ext}.bgRPGC.bs${bin_size}"
else
    bg_suffix=".bgRPGC.bs${bin_size}"
fi
bg_list=""
bg_out="${bg_prefix}${bg_suffix}"
for (( i=0; i<${#input_smooth_lengths[@]}; i++ )); do
    SMTH_LEN=${input_smooth_lengths[i]}
    WT=${input_smooth_weights[i]}
    bedgraph2="${bg_prefix}${bg_suffix}.smth${SMTH_LEN}.bedGraph"
    bg_list="${bg_list}${bedgraph2} "
    if [ ${#input_smooth_lengths[@]} -gt 1 ]; then
        if [ $i -eq 0 ]; then
            bg_out="${bg_out}.${WT}smth${SMTH_LEN}"
        else
            bg_out="${bg_out}_${WT}smth${SMTH_LEN}"
        fi
        outFormat="bedgraph"
        outExt=".bedGraph"
    else
        bg_out="${bg_out}.smth${SMTH_LEN}"
        outFormat="bigwig"
        outExt=".bw"
    fi
done


## Run bamCoverage on BAM file with the predetermined scale factor.

echo "Running deeptools bamCoverage on bam file..."

t0=$(date +%s)
bigwig1="${outdir}/$(basename --suffix='.bam' "${bam1}")${out_suffix}.bw"
bamCoverage -b "${bam1}" ${opts} \
    --scaleFactor ${RPGC_sf} \
    --extendReads ${read_ext} \
    --binSize ${bin_size} \
    --smoothLength ${smooth_length} \
    --normalizeUsing None \
    -p ${nproc} \
    -of bigwig \
    -o "${bigwig1}"
t1=$(date +%s)

echo "Finished bamCoverage in $((t1-t0)) s."

## Only generate smoothed Input tracks if they and the overall weight-averaged
## smoothed Input track (background/bias model) don't already exist.
if [ ! -f "${bg_out}.bw" ]; then
    for SMTH_LEN in ${input_smooth_lengths[@]} ; do
        ## Generate smoothed bedgraph/bigwig from bam if it doesn't already exist.
        bg2="${bg_prefix}${bg_suffix}.smth${SMTH_LEN}${outExt}"
        if [ ! -f "${bg2}" ]; then
            echo "Running bamCoverage using SMTH_LEN ${SMTH_LEN}..."
            t0=$(date +%s)
            bamCoverage -b "${bam1}" ${opts} \
                --scaleFactor ${RPGC_sf} \
                --extendReads ${read_ext} \
                --binSize ${bin_size} \
                --smoothLength ${SMTH_LEN} \
                --normalizeUsing None \
                -p ${nproc} \
                -of ${outFormat} \
                -o "${bg2}"
            t1=$(date +%s)
            echo "bamCoverage finished in $((t1-t0)) s."
            echo "outFormat is ${outFormat}."
            if [ "${outFormat}" = "bedgraph" ]; then
                echo "Sorting bedgraph using sort -k1,1 -k2,2n..."
                t0=$(date +%s)
                LC_ALL=C sort -k1,1 -k2,2n --parallel=8  -o "${bg2}" "${bg2}"
                t1=$(date +%s)
                echo "Finished sorting in $((t1-t0)) s."
            fi
        else
            echo "$(basename ${bg2}) already exists."
            echo "Skipping bamCoverage for SMTH_LEN ${SMTH_LEN}."
        fi
    done
else
    echo "Background/bias track $(basename ${bg_out}).bw already exists."
    echo "Using existing background/bias bigWig for ChIP-to-Input ratio."
fi

wait


## Construct a background/bias signal track by taking a weighted mean of
## the different smoothed pseudo-Input tracks (if it doesn't already exist).

if [ ! -f "${bg_out}.bw" ]; then

    echo "List of smoothed pseudo-background bedGraphs:"
    echo "   ${bg_list}"
    echo "Output name for weight-averaged smoothed pseudo-background bigWig:"
    echo "   $(basename ${bg_out}).bw"

    ## Compute weighted average of smoothed bg bedGraphs.
    echo "Computing (weighted) mean smoothed bg bigWig..."
    t0=$(date +%s)
    bedtools unionbedg -i ${bg_list} \
        | LC_ALL=C sort -k1,1 -k2,2n --parallel=8 \
        | awk -v W="${input_smooth_weights[*]}" 'BEGIN{split(W, weights); OFS="\t";}{
            s=0; w=0;
            for(i=1; i<=NF-3; i++){
                w+=weights[i];
                s+=weights[i]*$(i+3);
            }
            print $1, $2, $3, s/w; }' \
        > "${bg_out}.bedGraph"
    t1=$(date +%s)
    echo "Finished in $((t1-t0)) s."

    ## Convert the smoothed background/bias model from bedGraph to bigWig.
    echo "Converting $(basename ${bg_out}.bedGraph) to bigWig..."
    t0=$(date +%s)
    ${ucsc_dir}/bedGraphToBigWig "${bg_out}.bedGraph" "${chrom_sizes_path}" "${bg_out}.bw"
    t1=$(date +%s)
    echo "Finished in $((t1-t0)) s."

    ## Convert the individual smoothed pseduo-bg tracks to bigWig.
    ## Comment this block out if there is no need for these bigWigs.
    for f in ${bg_list}; do
        echo "Converting $f to bigWig..."
        t0=$(date +%s)
        bw_out="${outdir}/$(basename --suffix='.bedGraph' $f).bw"
        ${ucsc_dir}/bedGraphToBigWig "$f" "${chrom_sizes_path}" "${bw_out}"
        t1=$(date +%s)
        echo "Finished in $((t1-t0)) s."
    done

    ## Delete the unneeded bedGraph files.
    rm "${bg_out}.bedGraph"
    for f in ${bg_list}; do rm $f; done

else
    echo "Background/bias track $(basename ${bg_out}).bw already exists."
    echo "Using existing background/bias bigWig for ChIP-to-Input ratio."
fi

wait


## Run bigwigCompare to compute the ratio of normalized bam to smoothed bam pseudo-control.

echo 'Running bigwigCompare with "--operation ratio"...'
t0=$(date +%s)
ratioBigwig="${outdir}/$(basename --suffix='.bam' "${bam1}")"
if [ "${pseudocount}" != "0" ]; then
    ratioBigwig="${ratioBigwig}${out_suffix}.RPGC-SmthRPGC-Ratio_psdct${pseudocount}.bw"
else
    ratioBigwig="${ratioBigwig}${out_suffix}.ChIP-SmthChIP-Ratio.bw"
fi
bigwigCompare -b1 "${bigwig1}" -b2 "${bg_out}.bw" \
    --operation ratio \
    --pseudocount ${pseudocount} \
    --binSize ${bin_size} \
    -p ${nproc} \
    -of bigwig \
    -o "${ratioBigwig}"
t1=$(date +%s)
echo "Finished in $((t1-t0)) s."


## Also generate a bigWig file for the log2 ratio.

echo 'Running bigwigCompare with "--operation log2"...'
t0=$(date +%s)
log2Bigwig="${outdir}/$(basename --suffix='.bw' "${ratioBigwig}").log2.bw"
bigwigCompare -b1 "${bigwig1}" -b2 "${bg_out}.bw" \
    --operation log2 \
    --pseudocount ${pseudocount} \
    --binSize ${bin_size} \
    -p ${nproc} \
    -of bigwig \
    -o "${log2Bigwig}"
t1=$(date +%s)
echo "Finished in $((t1-t0)) s."

wait

echo "ALL DONE."


#############
###  END  ###
#############
