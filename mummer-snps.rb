#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2017-05-17
# version: 	0.01

require_relative 'genbank-manip'
# query_f=ARGV[0]
# ref_f=ARGV[1]
snps_file = ARGV[0]
gbk_file = ARGV[1]

def mummer_snps
  # nucmer --prefix=ref_qry ref.fasta qry.fasta
  # delta-filter -r -q ref_qry.delta > ref_qry.filter
  # show-snps -ClrT ref_qry.filter > ref_qry.snps
end


def read_genome_fts gbk

  fts_string = gbk.getFts

  features = []

  fts_string.split("\n").each do |l|
    lA = l.split("\t")
    lA[2] = lA[2].to_i
    lA[3] = lA[3].to_i
    features << lA
  end

  features

end


def build_gene_snps gbk, gbk_features, genes_mutations

  # [P1]    [SUB]   [SUB]   [P2]    [BUFF]  [DIST]  [LEN R] [LEN Q] [FRM]   [TAGS]
  # 4628    C       A       4628    4628    4628    6291459 6332858 1       1       PAE981  pae975-scaffold
  fout = File.open("zzz-output.tsv","w")
  header = [
    "Feature",
    "Strain",
    "Begin",
    "End",
    "Strand",
    "Partial",
    "Pseudo",
    "ProtID",
    "LocusTag",
    "Gene",
    "Prot-Product",
    "Nucleotide Mutation",
    "Protein Mutation"
  ]

  fout.write(header.join("\t"))
  fout.write("\n")

  seq = gbk.getSeq
  seq = seq.split("\n")[1..-1].join("")
  features = []

  genes_mutations.each do |ft_pos, mutation|

    next if ft_pos == "intergenic"

    features = []
    ft_pos_a = ft_pos.split("..")
    gene_seq_wt = seq[ft_pos_a[0].to_i-1..ft_pos_a[1].to_i-1]
    gene_seq_m = gene_seq_wt.clone

    nt_mutations = []

    mutation.reverse.each do |m|

      snp_gene_pos = m[0].to_i-ft_pos_a[0].to_i

      nt_mutations << "#{m[1]}#{snp_gene_pos}#{m[2]}"

      if m[1] == "."
        gene_seq_m.insert(snp_gene_pos, m[2])
      elsif m[2] == "."
        gene_seq_m.slice!(snp_gene_pos)
      else
        gene_seq_m[snp_gene_pos] = m[2]
      end

      if features.empty?
        features = m[12..-1]
      end

    end


    seq_aa_wt = Bio::Sequence.auto(gene_seq_wt)
    seq_aa_m = Bio::Sequence.auto(gene_seq_m)

    if ft_pos_a[2].to_i != 1
      seq_aa_wt.reverse_complement!
      seq_aa_m.reverse_complement!
    end

    seq_aa_wt = seq_aa_wt.translate
    seq_aa_m = seq_aa_m.translate

    bioalign = Bio::Alignment::OriginalAlignment.new([seq_aa_wt, seq_aa_m])
    mafft = Bio::MAFFT.new
    bioalign = bioalign.do_align(mafft)

    align_seqs = []
    bioalign.each_pair do |k, pairs|
      align_seqs << pairs
    end

    aa_mutations = []

    (align_seqs[0].length).times do |i|
      if align_seqs[0][i] != align_seqs[1][i]
        aa_mutations << align_seqs[0][i]+i.to_s+align_seqs[1][i]
        # puts "  AA change : " + align_seqs[0][i]+i.to_s+align_seqs[1][i]
      end
    end

    fout.write(features.join("\t"))
    fout.write("\t#{nt_mutations.join(',')}")
    fout.write("\t#{aa_mutations.join(',')}\n")
  end

  fout.close

end



def read_snps snps_file, gbk_features

  genes_mutations = {}
  genes_mutations["intergenic"] = []

  File.open(snps_file, "r") do |f|
    # skip header
    4.times{ |x| f.gets }
    while l=f.gets

      lA = l.chomp.split("\t")
      # puts "#{lA[0]}"
      snp_pos = lA[0].to_i

      mutation_information = []

      gbk_features.each do |ft|

        if ft[2] > snp_pos

          if mutation_information.empty?
            genes_mutations["intergenic"] << lA
          end

          break

        elsif snp_pos > ft[2] and snp_pos < ft[3]

          gene_key = "#{ft[2].to_i}..#{ft[3].to_i}..#{ft[4]}"
          if ! genes_mutations.has_key? gene_key
            genes_mutations[gene_key] = []
          end

          genes_mutations[gene_key] << lA + ft

        end

      end

      # puts lA.join("\t") + "\t" + mutation_information.join("\t")

    end
  end

  genes_mutations

end


# features = read_genome_fts genome_fts
gbk = GenbankParser.new(gbk_file)
gbk_features = read_genome_fts gbk

genes_mutations = read_snps snps_file, gbk_features
build_gene_snps gbk, gbk_features, genes_mutations
