#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	2017-11-01
# version: 	0.01

require 'parallel'
require_relative 'multifasta-manip'
require_relative 'genbank-manip'
require_relative 'cd-hit-parser'

# Need to have *.fasta files *.gb.fts.tsv

usage = "create-species-model.rb <dir> <threads>"


# Extract features from genbank files
def extract_gb_features file

  if ! File.exists? file+".fts.tsv"
    puts "  ..extracting features from #{gb_f}"
    gb_parser = GenbankParser.new(file)
    out = gb_parser.getFts
    File.open(file+".fts.tsv","w") do |fout|
      fout.write(out)
    end
  end

end


# Extract CDS, tRNA, rRNA and output sequences in corresponding files
# Skip partial genes
def extract_genes file_prefix

  if ! File.exists? file_prefix+".cds.fasta"

    puts "  ..extracting genes from #{file_prefix}"

    fasta_f = file_prefix + ".fasta"
    gbk_f = file_prefix + ".gb"
    gbk_parser = GenbankParser.new(gbk_f)
    cds = gbk_parser.getFtsProtSequences("#{File.basename(file_prefix)}", true)

    File.open(file_prefix+".cds.fasta","w") do |fout|
      fout.write(cds)
    end

    fts_f = file_prefix + ".gb.fts.tsv"
    fasta_parser = FastaParser.new(fasta_f)
    trna_output = ""
    rrna_output = ""

    File.open(fts_f, "r").drop(1).each do |f|
      lA = f.chomp!.split("\t")
      next if lA[5] == "1"        # skip partial gene
      basename = File.basename(fasta_f).gsub('.fasta','')
      header = "#{basename}|#{lA[1]}|#{lA[6]}|#{lA[7]}|#{lA[8]}|#{lA[9]}"
      if lA[0] == "tRNA"
        trna_output += fasta_parser.getSeqLoc(lA[1], "#{lA[2]}..#{lA[3]}", "#{lA[4]}", 0, header)
      elsif lA[0] == "rRNA"
        rrna_output += fasta_parser.getSeqLoc(lA[1], "#{lA[2]}..#{lA[3]}", "#{lA[4]}", 0, header)
      end
    end

    # write output
    File.open(file_prefix+".trna.fasta", "w") do |f_out|
      f_out.write(trna_output)
    end

    File.open(file_prefix+".rrna.fasta", "w") do |f_out|
      f_out.write(rrna_output)
    end

  end

end

# Parse cd_hit
def parse_cd_hit cdHitFile, fastaFile

  cdhit_parser = CdhitParser.new(cdHitFile, fastaFile)
  genome_index, genome_genes = cdhit_parser.build_genome_index
  clusters = cdhit_parser.parse_cdhit_file
  cdhit_parser.build_cluster_matrix(clusters, genome_index, genome_genes)
  cdhit_parser.build_cluster_info(clusters, genome_index)

end



def create_pangenome dir, proc

  Dir.mkdir("#{dir}/zz_pangenome") if ! File.exists? "#{dir}/zz_pangenome"

  files = [{in: "cds.fasta", out: "CDS.fasta"},
           {in: "trna.fasta", out: "tRNA.fasta"},
           {in: "rrna.fasta", out: "rRNA.fasta"}]

  Parallel.map(files, in_threads: proc) { |item|
    puts "  ..#{item[:in]}"
    `cat #{dir}/*.#{item[:in]} > #{dir}/zz_pangenome/#{item[:out]}`
    `cd-hit -c 0.80 -d 0 -M 0 -s 0.70 -i #{dir}/zz_pangenome/#{item[:out]} -o #{dir}/zz_pangenome/#{item[:out]}.cd-hit`
    parse_cd_hit("#{dir}/zz_pangenome/#{item[:out]}.cd-hit.clstr", "#{dir}/zz_pangenome/#{item[:out]}")
  }

end



# Main #
dir = ARGV[0]
proc = ARGV[1]
proc = proc.to_i

# Extract Features from Genbank
all_gb = Dir.glob("#{dir}/*.gb")
all_gb.delete("#{dir}/OTHER_SAMPLES.gb")
puts "## Extracting Features "
Parallel.map(all_gb, in_processes: proc) { |gb_f|
  extract_gb_features gb_f
}

# Extract Genes from Sample Features
all_sample_fts = Dir.glob("#{dir}/*.gb.fts.tsv")
puts "## Extracting Genes "
Parallel.map(all_sample_fts, in_processes: proc) { |fts_f|
  file_prefix = fts_f.gsub(".gb.fts.tsv","")
  extract_genes file_prefix
}

# Creating Pan-Genome
puts "## Creating Pan-Genome model "
create_pangenome dir, proc
