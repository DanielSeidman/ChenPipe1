#expects $1 to be a list of sample directories
#expects $2 to be the name of the sample file you want to make


printf "BioSample,Run,fq1,fq2\n" > ${2};
for pathval in $(cat ${1});
  do ls $pathval/*1.fq.gz | awk 'BEGIN{
                                  OFS=","
                                }{
                                  dir_n=split($1,directories,"/");
                                  token_n=split(directories[dir_n],tokens,"_");
                                  split($1,no_suffix,"1.fq.gz");
                                  print tokens[1],tokens[1]"_"tokens[2]"_"tokens[3]"_"tokens[4],no_suffix[1]"1.fq.gz",no_suffix[1]"2.fq.gz"
                                }' >> ${2};
  done
