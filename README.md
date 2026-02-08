**Required Python libraries**

Python libraries: `numpy`, `scipy`, `matplotlib`, `pandas`.

**TM-align and DALI alignments for SCOP40**

These alignments should be installed from the Foldseek supp. data.

<pre>
wget https://wwwuser.gwdg.de/~compbiol/foldseek/scop40pdb.tar.gz
tar -zxf scop40pdb.tar.gz
cp alignResults/rawoutput/tmaln $MY_REPOS/statsig_paper/big_data/tm.scop40
cp alignResults/rawoutput/dalialn $MY_REPOS/statsig_paper/big_data/dalis.scop40
</pre>

Create Foldseek SCOP40 database in `../big_foldseek/scop40/scop40db`.
