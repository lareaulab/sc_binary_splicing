# Running rMATS on Trapnell
# It is unstranded, thus no strand specificity

DIR="/mnt/lareaulab/cfbuenabadn/RNASeq/Human/Trapnell/SRP033135" 
OUT="/mnt/lareaulab/cfbuenabadn/RNASeq/Human/Trapnell/splicing/rMATS/output" 
BAM="star_output/Aligned.sortedByCoord.out.bam"
mkdir output

# H00 vs H24
mkdir $OUT/H24
python ~/software/rMATS.3.2.5/RNASeq-MATS.py -b1 $DIR/SRR1033282_1/$BAM,$DIR/SRR1033283_1/$BAM,$DIR/SRR1033284_1/$BAM \
                                             -b2 $DIR/SRR1033285_1/$BAM,$DIR/SRR1033286_1/$BAM,$DIR/SRR1033287_1/$BAM \
                                             -gtf ~/Genomes/Human/hg38/gencode.v24.annotation.gtf \
                                             -o $OUT/H24 \
                                             -t paired \
                                             -len 100 \
                                             -c 0.0001

# H00 vs H48
mkdir $OUT/H48
python ~/software/rMATS.3.2.5/RNASeq-MATS.py -b1 $DIR/SRR1033282_1/$BAM,$DIR/SRR1033283_1/$BAM,$DIR/SRR1033284_1/$BAM \
                                             -b2 $DIR/SRR1033288_1/$BAM,$DIR/SRR1033289_1/$BAM,$DIR/SRR1033290_1/$BAM \
                                             -gtf ~/Genomes/Human/hg38/gencode.v24.annotation.gtf \
                                             -o $OUT/H48 \
                                             -t paired \
                                             -len 100 \
                                             -c 0.0001

# H00 vs H72
mkdir $OUT/H72
python ~/software/rMATS.3.2.5/RNASeq-MATS.py -b1 $DIR/SRR1033282_1/$BAM,$DIR/SRR1033283_1/$BAM,$DIR/SRR1033284_1/$BAM \
                                             -b2 $DIR/SRR1033291_1/$BAM,$DIR/SRR1033292_1/$BAM,$DIR/SRR1033293_1/$BAM \
                                             -gtf ~/Genomes/Human/hg38/gencode.v24.annotation.gtf \
                                             -o $OUT/H72 \
                                             -t paired \
                                             -len 100 \
                                             -c 0.0001



