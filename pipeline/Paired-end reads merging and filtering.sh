## 1.Remove barcodes and primers.  
    cd /mnt/...
    
    mkdir fastp

    ## V3V4
    for filename in *.1.fq.gz
    do
    base=$(basename $filename .1.fq.gz)
    echo $base

    fastp -i ${base}.1.fq.gz -I ${base}.2.fq.gz \
      -o fastp/${base}.1.fastp.fq.gz -O fastp/${base}.2.fastp.fq.gz \
      -w 6 \
      -A \
      -f 26\
      -F 26 \
      -t 25 \
      -T 25 \
      -l 150 \
      -n 0\
      -Q
    done

    ## V4

    mkdir fastp

    for filename in *_1.fastq.gz
    do
    base=$(basename $filename _1.fastq.gz)
    echo $base

    fastp -i ${base}_1.fastq.gz -I ${base}_2.fastq.gz \
      -o fastp/${base}.1.fastp.fq.gz -O fastp/${base}.2.fastp.fq.gz\
      -w 6 \
      -A \
      -f 27\
      -F 28 \
      -t 23 \
      -T 22 \
      -l 100 \
      -n 0\
      -Q
    done
    
## 2. Merge paired-end reads.

    # V3V4
    mkdir flash

    for filename in *.1.fastp.fq.gz
    do
    base=$(basename $filename .1.fastp.fq.gz)
    echo $base

    flash \
     ${base}.1.fastp.fq.gz ${base}.2.fastp.fq.gz \
    -p 33 -x 0.25 -m 30  -M 193 -t 6 \
    -o flash/${base}
    done

    # V4
    mkdir flash

    for filename in *.1.fastp.fq.gz
    do
    base=$(basename $filename .1.fastp.fq.gz)
    echo $base

    flash \
     ${base}.1.fastp.fq.gz ${base}.2.fastp.fq.gz \
    -p 33 -x 0.25 -m 30  -M 199 -t 6 \
    -o flash/${base}
    done

   
## 3. Filter sequence by length.
   cd /mnt/...

    # V3V4
    mkdir flash_qc

    for filename in *.flash.fq.gz
    do
    base=$(basename $filename .flash.fq.gz)
    echo $base

    fastp -i ${base}.flash.fq.gz \
      -o flash_qc/${base}.qc.flash.fq.gz \
      -w 6 \
      -A \
      -l 400\
      -Q
    done
    
    #V4
    mkdir flash_qc

    for filename in *.flash.fq.gz
    do
    base=$(basename $filename .flash.fq.gz)
    echo $base

    fastp -i ${base}.flash.fq.gz \
      -o flash_qc/${base}.qc.flash.fq.gz \
      -w 6 \
      -A \
      -l 200\
      -Q
    done
