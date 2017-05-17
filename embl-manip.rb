#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# genbank-manip.rb --- 
# 
# Filename: embl-manip.rb
# Description: 
# Author: Maxime Déraspe
# Maintainer: 
# Created: Tue Nov  6 10:40:15 2012 (-0500)
# Version: 
# Last-Updated: 
#           By: 
#     Update #: 0
# Fonctions available TODO :
#
# 	JUST A COPY OF genbank-manip.rb for embl
# Code:
require 'bio'


# Help Usage
usage = 
"Usage : $ embl-manip.rb <option> <parameters> <gbk file>

	options:
		getfts 	return all fts info
		getfts-dna	return all fts dna sequences
		getfts-pep	return all fts pep sequences
		getft  		<kword> [pep] nt by default
		getft-loc	beg..end

"

## Init ##

opt = ARGV[0].to_s.downcase
abort "File Not Found ! \n\n#{usage}" if ! File.exist?(ARGV[ARGV.length-1].to_s)
f = ARGV[ARGV.length-1] or abort usage
@gbfile = File.new(f, "r").read
@gb = Bio::EMBL.new(@gbfile)
@accession = @gb.accession
@org = @gb.species

## Functions ##

# Fct: Get dna sequence
def getDna (cds, seq)
  loc = cds.locations
  sbeg = loc[0].from.to_i
  send = loc[0].to.to_i
  fasta = Bio::Sequence::NA.new(seq.subseq(sbeg,send))
  position = "#{sbeg}..#{send}"
  if loc[0].strand == -1
    fasta.reverse_complement!
  end
  dna = Bio::Sequence.auto(fasta)
  return dna
end

# Fct: Get peptides sequence
def getPep (cds, seq)

end

# Fct: Get Features Sequences
def getFtsSequences
  @gb.each_cds do |ft|
    ftH = ft.to_hash
    loc = ft.locations
    loc = "c#{ft.locations[0].to_s}" if ft.locations[0].strand == -1
    gene = []
    product = []
    gene = ftH["gene"] if !ftH["gene"].nil?
    product = ftH["product"] if !ftH["product"].nil?
    dna = getDna(ft,@gb.to_biosequence)
    seqout = dna.output_fasta("#{@accession}|#{loc}|#{ftH["protein_id"][0]}|#{gene[0]}|#{product[0]}|#{@org}",60)
    puts seqout
  end
end

# Fct: Get all features location
def getFtsLoc
  location = ARGV[1]
  loc = location.split("..")
  @gb.each_cds do |ft|
    ftH = ft.to_hash
    ftloc = ft.locations
    if ftloc[0].from == loc[0].to_i && ftloc[0].to == loc[1].to_i
      gene = []
      product = []
      gene = ftH["gene"] if !ftH["gene"].nil?
      product = ftH["product"] if !ftH["product"].nil?
      loc = "c#{location}" if ftloc[0].strand == -1
      dna = getDna(ft,@gb.to_biosequence)
      seqout = dna.output_fasta("#{@accession}|#{loc}|#{ftH["protein_id"][0]}|#{gene[0]}|#{product[0]}|#{@org}",60)
      puts seqout
    end
  end
end

# Fct: Get feature of according gene or product matching
#      exactly the kword (insensitive)
def getFt 
  kword = ARGV[1]
  seq = @gb.to_biosequence
  seqoptions = ""

  for c in 2..ARGV.length-1
    seqoptions += "#{ARGV[c]},"
  end
  
  # look through all features
  @gb.each_cds do |ft|
    ftH = ft.to_hash
    loc = ft.locations
    gene = []
    product = []
    if (!ftH["gene"].nil? && ftH["gene"][0].downcase.include?(kword.downcase)) or
        (!ftH["product"].nil? && ftH["product"][0].downcase.include?(kword.downcase))      
      sbeg = loc[0].from.to_i
      send = loc[0].to.to_i
      fasta = Bio::Sequence::NA.new(seq.subseq(sbeg,send))
      position = "#{sbeg}..#{send}"
      if loc[0].strand == -1
        fasta.reverse_complement!
        position = "c#{position}"
      end
      pep = Bio::Sequence.new(fasta.translate)
      gene = ftH["gene"][0] if !ftH["gene"].nil?
      product = ftH["product"][0] if !ftH["product"].nil?
      if seqoptions.downcase.include?("pep") or seqoptions.downcase.include?("prot")
        puts pep.output_fasta("#{@accession}|#{position}|#{ftH["protein_id"][0]}|#{gene}|#{product}|#{@org}", 60)
      else
        dna = Bio::Sequence.auto(fasta)
        puts dna.output_fasta("#{@accession}|#{position}|#{ftH["protein_id"][0]}|#{gene}|#{product}|#{@org}",60)
      end
    end
  end
end


# Fct: Get all features info
def getFts
  puts "CDS\tLocation\tGene\tProduct"
  @gb.each_cds do |ft|
    ftH = ft.to_hash
    loc = ft.locations
    seqBeg = loc[0].from.to_s
    seqEnd = loc[0].to.to_s
    strand = loc[0].strand.to_s
    gene = []
    product = []
    gene = ftH["gene"] if !ftH["gene"].nil?
    product = ftH["product"] if !ftH["product"].nil?
    protId = ftH["protein_id"][0] if !ftH["protein_id"].nil?
    puts "#{@accession}\t#{seqBeg}\t#{seqEnd}\t#{strand}\t#{protId}\t#{gene[0]}\t#{product[0]}"
  end
end


# Fct: Merge genbank annotation to fasta Sequence
def mergeGbkSeq
  
end

# Fct: Shift genbank annotation to fasta sequence
def shiftGbk

end


## Main ##

# switch case with the options
case opt

when "getfts"
  getFts
when "getfts-dna"
  getFtsSequences
when "getft-loc"
  getFtsLoc
when "getft"
  getFt
else
  puts usage
end

