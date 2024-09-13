#used to count allele ratios for a single individual
#useful for checking the ploidy/contamination
#samtools mpileup -q 10 -Q 10 -s bam1 [bam2 ...]|awk -f pileupCounts.awk
BEGIN{
	if (limit3 == "")
		limit3 = 5 #minimum allele coverage to consider variant
	if (limit4 == "")
		limit4 = 30 #minimum total coverage
	if (limit5 == "")
		limit5 =100 #maximum total coverage
	FS="\t"
	OFS="\t"
}

function abs(a)
{
	if (a < 0)
		return -a
	else
		return a
}

{
	sum = 0
	for (i = 4; i <= NF; i+=4)
		sum += $i;

	if (sum < limit4 || sum > limit5)
		next

	delete c
	for (i = 5; i <= NF; i+=4) {
		if ($(i-1) == 0)
			$i = ""
		gsub(/\$/,"",$i)  #remove end of reads
		gsub(/\^./,"",$i) #remove quality
		while (match($i, /[+-][1-9][0-9]*/) > 0) { #remove indels
			$i = substr($i, 1, RSTART - 1) substr($i, RSTART + RLENGTH + substr($i, RSTART + 1, RLENGTH - 1))
		}
		tmp = $i
		c[0] += gsub(/[Aa]/, "", tmp)
		c[1] += gsub(/[Cc]/, "", tmp)
		c[2] += gsub(/[Gg]/, "", tmp)
		c[3] += gsub(/[Tt]/, "", tmp)
	}
	alleles = 0
	sum = 0
	if (c[0] >= limit3) {
		++alleles
		sum+=c[0]
	}
	if (c[1] >= limit3) {
		++alleles
		sum+=c[1]
	}
	if (c[2] >= limit3) {
		++alleles;
		sum+=c[2]
	}
	if (c[3] >= limit3) {
		sum+=c[3]
		++alleles;
	}

#	if (alleles >= 2)
#		print c[0] "\t" c[1] "\t" c[2] "\t" c[3]
	if (alleles == 2 && sum >= limit4) {
		printf $1 "\t" $2 "\t"
		for (i = 0; i < 4; ++i)
			if (c[i] >= limit3)
				printf c[i] "\t"
		print ""
	}
		
}
