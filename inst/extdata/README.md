#Liftover EAS_fourier_ls-all.bed to hg38
Wrote lifted bedfile file to XPASS/inst/extdata/EAS_fourier_ls-all_Hg38.bed. 1442 out of 1445 were lifted.

#LiftOver EUR_fourier_ls-all.bed to hg38
Wrote lifted bedfile file to XPASS/inst/extdata/EUR_fourier_ls-all_Hg38.bed. 1700 out of 1703 were lifted.

wget https://bitbucket.org/nygcresearch/ldetect-data/get/ac125e47bf7f.zip
unzip ac125e47bf7f.zip
rm ac125e47bf7f.zip
cp nygcresearch-ldetect-data-ac125e47bf7f/AFR/fourier_ls-all.bed AFR_fourier_ls-all.bed
rm -rf nygcresearch-ldetect-data-ac125e47bf7f/

##LiftOver EUR_fourier_ls-all.bed to hg38
Wrote lifted bedfile file to XPASS/inst/extdata/AFR_fourier_ls-all_Hg38.bed. 2580 out of 2581 were lifted.


#in R

```
#liftover(x="/home/bbita/R/x86_64-pc-linux-gnu-library/4.0/XPASS/extdata/EAS_fourier_ls-all.bed", pop='EAS', hg_old="hg19", hg_new="hg38")
#liftover(x="/home/bbita/R/x86_64-pc-linux-gnu-library/4.0/XPASS/extdata/EUR_fourier_ls-all.bed", pop='EUR', hg_old="hg19", hg_new="hg38")
#liftover(x="XPASS/inst/extdata/AFR_fourier_ls-all.bed", pop='AFR', hg_old="hg19", hg_new="hg38")
```

