# Running rMATS on Hubbard

DIR="/mnt/lareaulab/cfbuenabadn/RNASeq/Mouse/HubbardRNASeq/SRP017778" 
OUT="/mnt/lareaulab/cfbuenabadn/RNASeq/Mouse/HubbardRNASeq/splicing/rMATS/output" 
BAM="star_output/Aligned.sortedByCoord.out.bam"
#mkdir output

# -8D vs -4D
mkdir $OUT/ND4
python ~/software/rMATS.3.2.5/RNASeq-MATS.py -b1 $DIR/SRR645824_1/$BAM,$DIR/SRR645826_1/$BAM,$DIR/SRR645828_1/$BAM,$DIR/SRR645830_1/$BAM \
                                             -b2 $DIR/SRR645833_1/$BAM,$DIR/SRR645835_1/$BAM,$DIR/SRR645837_1/$BAM \
                                             -gtf ~/Genomes/Mouse/Gencode/gencode.vM10.annotation.gtf \
                                             -o $OUT/DN4 \
                                             -t paired \
                                             -len 50 \
                                             -c 0.0001

# -8D vs D0
mkdir $OUT/D0
python ~/software/rMATS.3.2.5/RNASeq-MATS.py -b1 $DIR/SRR645824_1/$BAM,$DIR/SRR645826_1/$BAM,$DIR/SRR645828_1/$BAM,$DIR/SRR645830_1/$BAM \
                                             -b2 $DIR/SRR645840_1/$BAM,$DIR/SRR645842_1/$BAM,$DIR/SRR645844_1/$BAM \
                                             -gtf ~/Genomes/Mouse/Gencode/gencode.vM10.annotation.gtf \
                                             -o $OUT/D0 \
                                             -t paired \
                                             -len 50 \
                                             -c 0.0001

# -8D vs D1
mkdir $OUT/D1
python ~/software/rMATS.3.2.5/RNASeq-MATS.py -b1 $DIR/SRR645824_1/$BAM,$DIR/SRR645826_1/$BAM,$DIR/SRR645828_1/$BAM,$DIR/SRR645830_1/$BAM \
                                             -b2 $DIR/SRR645930_1/$BAM,$DIR/SRR645932_1/$BAM,$DIR/SRR645934_1/$BAM,$DIR/SRR646092_1/$BAM \
                                             -gtf ~/Genomes/Mouse/Gencode/gencode.vM10.annotation.gtf \
                                             -o $OUT/D1 \
                                             -t paired \
                                             -len 50 \
                                             -c 0.0001

# -8D vs D7
mkdir $OUT/D7
python ~/software/rMATS.3.2.5/RNASeq-MATS.py -b1 $DIR/SRR645824_1/$BAM,$DIR/SRR645826_1/$BAM,$DIR/SRR645828_1/$BAM,$DIR/SRR645830_1/$BAM \
                                             -b2 $DIR/SRR645846_1/$BAM,$DIR/SRR645849_1/$BAM,$DIR/SRR645851_1/$BAM,$DIR/SRR645853_1/$BAM,$DIR/SRR646182_1/$BAM \
                                             -gtf ~/Genomes/Mouse/Gencode/gencode.vM10.annotation.gtf \
                                             -o $OUT/D7 \
                                             -t paired \
                                             -len 50 \
                                             -c 0.0001

# -8D vs D16
mkdir $OUT/D16
python ~/software/rMATS.3.2.5/RNASeq-MATS.py -b1 $DIR/SRR645824_1/$BAM,$DIR/SRR645826_1/$BAM,$DIR/SRR645828_1/$BAM,$DIR/SRR645830_1/$BAM \
                                             -b2 $DIR/SRR645855_1/$BAM,$DIR/SRR645857_1/$BAM,$DIR/SRR645859_1/$BAM,$DIR/SRR645861_1/$BAM \
                                             -gtf ~/Genomes/Mouse/Gencode/gencode.vM10.annotation.gtf \
                                             -o $OUT/D16 \
                                             -t paired \
                                             -len 50 \
                                             -c 0.0001

# -8D vs D21
mkdir $OUT/D21
python ~/software/rMATS.3.2.5/RNASeq-MATS.py -b1 $DIR/SRR645824_1/$BAM,$DIR/SRR645826_1/$BAM,$DIR/SRR645828_1/$BAM,$DIR/SRR645830_1/$BAM \
                                             -b2 $DIR/SRR645863_1/$BAM,$DIR/SRR645865_1/$BAM,$DIR/SRR645867_1/$BAM,$DIR/SRR645870_1/$BAM \
                                             -gtf ~/Genomes/Mouse/Gencode/gencode.vM10.annotation.gtf \
                                             -o $OUT/D21 \
                                             -t paired \
                                             -len 50 \
                                             -c 0.0001

# -8D vs D28
mkdir $OUT/D28
python ~/software/rMATS.3.2.5/RNASeq-MATS.py -b1 $DIR/SRR645824_1/$BAM,$DIR/SRR645826_1/$BAM,$DIR/SRR645828_1/$BAM,$DIR/SRR645830_1/$BAM \
                                             -b2 $DIR/SRR645872_1/$BAM,$DIR/SRR645875_1/$BAM,$DIR/SRR645877_1/$BAM,$DIR/SRR645879_1/$BAM \
                                             -gtf ~/Genomes/Mouse/Gencode/gencode.vM10.annotation.gtf \
                                             -o $OUT/D28 \
                                             -t paired \
                                             -len 50 \
                                             -c 0.0001

