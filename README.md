# Seqpare

Seqpare is a self-consistent metric and tool for comparing sequences based on the total effective overlaps of their interval sets. With this metric, the similarity of two interval sets is quantified by a single index, which directly represents the percentage of their effective overlapping: a similarity index of zero indicates totally unrelated sequences, and an index of one means that the interval sets are exactly the same. Seqpare tool provides functions for both searching and mapping large-scale interval datasets. The seqpare code is based on [ailist][ailist].

## Installing Seqpare

```
git clone https://github.com/deepstanding/seqpare
cd seqpare
make
sudo cp bin/seqpare /usr/local/bin
```
## Running Seqpare

### Compare two interval sets

```
seqpare <BED file 1> <BED file 2> -o output
```

### Search: compare a query interval set over a collection of interval sets

```
seqpare <BED file folder> <query BED file> -m 1 -o output
```

### Similarity mapping of two different collections of interval sets  

```
seqpare <BED file folder> <query BED file> -m 2 -o output
```

### Self similarity mapping of a collection of interval sets

```
seqpare <BED file folder> -m 3 -o output
```

## Data source

http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/

## Test results

Four files in the Test_results folder.

[ailist]: https://github.com/databio/AIList
