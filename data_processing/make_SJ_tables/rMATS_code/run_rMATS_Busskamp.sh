# Running rMATS on Velasco

DIR="/mnt/lareaulab/cfbuenabadn/RNASeq/Human/Busskamp/SRP045632" 
OUT="/mnt/lareaulab/cfbuenabadn/RNASeq/Human/Busskamp/splicing/rMATS/output" 
BAM="star_output/Aligned.sortedByCoord.out.bam"
mkdir output

# D0 vs D1
mkdir $OUT/D1
python ~/software/rMATS.3.2.5/RNASeq-MATS.py -b1 $DIR/SRR1554436_1/$BAM,$DIR/SRR1554437_1/$BAM,$DIR/SRR1554438_1/$BAM \
                                             -b2 $DIR/SRR1554439_1/$BAM,$DIR/SRR1554440_1/$BAM,$DIR/SRR1554441_1/$BAM \
                                             -gtf ~/Genomes/Human/gencode.v24.annotation.gtf \
                                             -o $OUT/D1 \
                                             -t paired \
                                             -len 101 \
                                             -c 0.0001

# D0 vs D3
mkdir $OUT/D3
python ~/software/rMATS.3.2.5/RNASeq-MATS.py -b1 $DIR/SRR1554436_1/$BAM,$DIR/SRR1554437_1/$BAM,$DIR/SRR1554438_1/$BAM \
                                             -b2 $DIR/SRR1554442_1/$BAM,$DIR/SRR1554443_1/$BAM,$DIR/SRR1554444_1/$BAM \
                                             -gtf ~/Genomes/Human/gencode.v24.annotation.gtf \
                                             -o $OUT/D3 \
                                             -t paired \
                                             -len 101 \
                                             -c 0.0001

# D0 vs D4
mkdir $OUT/D4
python ~/software/rMATS.3.2.5/RNASeq-MATS.py -b1 $DIR/SRR1554436_1/$BAM,$DIR/SRR1554437_1/$BAM,$DIR/SRR1554438_1/$BAM \
                                             -b2 $DIR/SRR1554445_1/$BAM,$DIR/SRR1554446_1/$BAM,$DIR/SRR1554447_1/$BAM \
                                             -gtf ~/Genomes/Human/gencode.v24.annotation.gtf \
                                             -o $OUT/D4 \
                                             -t paired \
                                             -len 101 \
                                             -c 0.0001



