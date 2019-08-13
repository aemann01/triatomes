#filter out tree that you want with archaeoptryx
#pull unaligned reads from LTP database
cat query.ids | while read line; do grep -w $line LTPs132_datasets.fasta/LTPs132_SSU_compressed.fasta -A 1 ; done > LTPs132_salmonella_ref.fa

#change uracil to thymine
sed -i '/^[^>]/ y/uU/tT/' LTPs132_datasets.fasta/LTPs132_SSU_aligned.fasta

#remove word wrap
remove_wordwrap.sh LTPs132_datasets.fasta/LTPs132_SSU_aligned.fasta > LTPs132_datasets.fasta/temp
mv LTPs132_datasets.fasta/temp LTPs132_datasets.fasta/LTPs132_SSU_aligned.fasta

#filter reference alignment
cat query.ids | while read line; do grep -w $line LTPs132_datasets.fasta/LTPs132_SSU_aligned.fasta -A 1 ; done > LTPs132_salmonella_ref.align.fa

#build tree
align_seqs.py -i ASV18.fa -t LTPs132_salmonella_ref.align.fa -o htes_align -p 80
cat htes_align/ASV18_aligned.fasta LTPs132_salmonella_ref.align.fa > query_ref.align.fa
sed -i 's/\./-/g' query_ref.align.fa
raxmlHPC-PTHREADS-SSE3 -f a -N 100 -G 0.2 -m GTRCAT -n constrain.tre -s queryalign_plus_refalign.fa -g pruned_tree.tre -T 4 -x 25734 -p 78543