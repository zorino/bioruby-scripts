#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# genbank-manip.rb --- 
# 
# Filename: genbank-manip.rb
# Description: 
# Author: Maxime Déraspe
# Maintainer: 
# Created: Tue Nov  6 10:40:15 2012 (-0500)
# Version: 
# Last-Updated: 
#           By: 
#     Update #: 0
#
# Fonctions TODO :
#
# Shift genbank file (by X nucleotides)
# Change orientation (strand) of genbank file
# Merge genbank files into one
#    Being able to combine those three functions is a must

# Code:

require 'bio'


## Init ##


class GenbankManipulator

  attr_accessor :gbkObj, :accession, :org

  def initialize (gbFile)
    @gbkObj = Bio::GenBank.new(gbFile)
    @accession = @gbkObj.accession
    @org = @gbkObj.organism
  end

  #Fct: Get Dna Sequence
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

  # Fct: Get Features Sequences
  def getFtsNtSequences
    @gbkObj.each_cds do |ft|
      ftH = ft.to_hash
      loc = ft.locations
      gene = []
      product = []
      protId = ""
      gene = ftH["gene"] if !ftH["gene"].nil?
      product = ftH["product"] if !ftH["product"].nil?
      protId = ftH["protein_id"][0] if !ftH["protein_id"].nil?
      locustag = ftH["locus_tag"][0] if !ftH["locus_tag"].nil?
      dna = getDna(ft,@gbkObj.to_biosequence)
      seqout = dna.output_fasta("#{@accession}|#{loc}|#{protId}|#{locustag}|#{gene[0]}|#{product[0]}",60)
      puts seqout
    end
  end

  # Fct: Get Features Sequences
  def getFtsProtSequences
    @gbkObj.each_cds do |ft|
      ftH = ft.to_hash
      loc = ft.locations
      gene = []
      product = []
      protId = ""
      if ftH.has_key? "pseudo"
        next
      end
      gene = ftH["gene"] if !ftH["gene"].nil?
      product = ftH["product"] if !ftH["product"].nil?
      protId = ftH["protein_id"][0] if !ftH["protein_id"].nil?
      locustag = ftH["locus_tag"][0] if !ftH["locus_tag"].nil?
      dna = getDna(ft,@gbkObj.to_biosequence)
      pep = ftH["translation"][0] if !ftH["translation"].nil?
      pepBioSeq = Bio::Sequence.auto(pep)
      seqout = pepBioSeq.output_fasta("#{@accession}|#{loc}|#{protId}|#{locustag}|#{gene[0]}|#{product[0]}",60)
      puts seqout
    end
  end


  # Fct: Get all features location
  def getFtsLoc
    location = ARGV[1]
    loc = location.split("..")
    protId = ""
    @gbkObj.each_cds do |ft|
      ftH = ft.to_hash
      ftloc = ft.locations
      if ftloc[0].from == loc[0].to_i && ftloc[0].to == loc[1].to_i
        gene = []
        product = []
        gene = ftH["gene"] if !ftH["gene"].nil?
        product = ftH["product"] if !ftH["product"].nil?
        protId = ftH["protein_id"][0] if !ftH["protein_id"].nil?
        location = "c#{location}" if ftloc[0].strand == -1
        dna = getDna(ft,@gbkObj.to_biosequence)
        seqout = dna.output_fasta("#{@accession}|#{location}|#{protId}|#{gene[0]}|#{product[0]}",60)
        puts seqout
      end
    end
  end

  # Fct: Get a feature at a particular locus_tag
  def getFtLocus
    locustag = ARGV[1]
    protId = ""
    @gbkObj.each_cds do |ft|
      ftH = ft.to_hash
      ftloc = ft.locations
      if ftH["locus_tag"][0].eql? locustag
        gene = []
        product = []
        gene = ftH["gene"] if !ftH["gene"].nil?
        product = ftH["product"] if !ftH["product"].nil?
        protId = ftH["protein_id"][0] if !ftH["protein_id"].nil?
        locustag = ftH["locus_tag"][0] if !ftH["locus_tag"].nil?
        if ftloc[0].strand == -1
          location = "c#{ftloc[0].from}..#{ftloc[0].to}"
        else
          location = "#{ftloc[0].from}..#{ftloc[0].to}"
        end
        dna = getDna(ft,@gbkObj.to_biosequence)
        seqout = dna.output_fasta("#{@accession}|#{location}|#{protId}|#{locustag}|#{gene[0]}|#{product[0]}",60)
        puts seqout
      end
    end
  end

  # Fct: get ft dna for a particular protein_id
  def getFtProtID
    protein_id = ARGV[1]
    protId = ""
    @gbkObj.each_cds do |ft|
      ftH = ft.to_hash
      ftloc = ft.locations
      if ftH["protein_id"][0].include? protein_id
        gene = []
        product = []
        gene = ftH["gene"] if !ftH["gene"].nil?
        product = ftH["product"] if !ftH["product"].nil?
        protId = ftH["protein_id"][0] if !ftH["protein_id"].nil?
        locustag = ftH["locus_tag"][0] if !ftH["locus_tag"].nil?
        if ftloc[0].strand == -1
          location = "c#{ftloc[0].from}..#{ftloc[0].to}"
        else
          location = "#{ftloc[0].from}..#{ftloc[0].to}"
        end
        dna = getDna(ft,@gbkObj.to_biosequence)
        seqout = dna.output_fasta("#{@accession}|#{location}|#{protId}|#{locustag}|#{gene[0]}|#{product[0]}",60)
        puts seqout
      end
    end
  end


  # Fct: Return the full dna sequence
  def getSeq
    bioSeq = @gbkObj.to_biosequence
    sequence = Bio::Sequence.new(bioSeq)
    puts sequence.output_fasta("#{bioSeq.accessions[0]}",60)
  end


  # Fct: Return part of the sequence specify by loc
  def getSeqLoc
    if ARGV.length == 4
      strand = ARGV[2]
      location = ARGV[1]
    else
      abort "You need to specify location and strand !"
    end
    loc = location.split("..")
    bioSeq = @gbkObj.to_biosequence
    if strand.to_i == -1
      sequence = Bio::Sequence.new(bioSeq.subseq(loc[0].to_i,loc[1].to_i).reverse_complement)
    elsif strand.to_i == 1
      sequence = Bio::Sequence.new(bioSeq.subseq(loc[0].to_i,loc[1].to_i))
    else
      abort "Bad Strand : 1 or -1 needed"
    end
    puts sequence.output_fasta("#{bioSeq.accessions[0]}|#{loc[0]}..#{loc[1]}|#{strand}",60)
  end


  # Fct: Get feature of according gene or product matching
  #      exactly the kword (insensitive)
  def getFt 
    kword = ARGV[1]
    seq = @gbkObj.to_biosequence
    seqoptions = ""  
    for c in 2..ARGV.length-1
      seqoptions += "#{ARGV[c]},"
    end
    # look through all features
    @gbkObj.each_cds do |ft|
      ftH = ft.to_hash
      loc = ft.locations
      gene = []
      product = []
      protId = ""
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
        protId = ftH["protein_id"][0] if !ftH["protein_id"].nil?
        locustag = ftH["locus_tag"][0] if !ftH["locus_tag"].nil?
        if seqoptions.downcase.include?("pep") or seqoptions.downcase.include?("prot")
          puts pep.output_fasta("#{@accession}|#{position}|#{protId}|#{locustag}|#{gene}|#{product}", 60)
        else
          dna = Bio::Sequence.auto(fasta)
          puts dna.output_fasta("#{@accession}|#{position}|#{protId}|#{locustag}|#{gene}|#{product}",60)
        end
      end
    end
  end


  # Fct: Get all features info
  def getFts
    puts "CDS\tLocation\tGene\tProduct"
    @gbkObj.each_cds do |ft|
      ftH = ft.to_hash
      loc = ft.locations
      seqBeg = loc[0].from.to_s
      seqEnd = loc[0].to.to_s
      strand = loc[0].strand.to_s
      gene = []
      product = []
      protId = ""
      gene = ftH["gene"] if !ftH["gene"].nil?
      product = ftH["product"] if !ftH["product"].nil?
      protId = ftH["protein_id"][0] if !ftH["protein_id"].nil?
      locustag = ftH["locus_tag"][0] if !ftH["locus_tag"].nil?
      puts "#{@accession}\t#{seqBeg}\t#{seqEnd}\t#{strand}\t#{protId}\t#{locustag}\t#{gene[0]}\t#{product[0]}"
    end
  end

  # Fct: GetTaxon out of GenbankFile
  def getTaxon
    classification = ['Kingdom:', 'Phylum:', 'Class:', 'Order:', 'Family:', 'Genus:', 'Species:', 'Organism:']
    i = 0
    @gbkObj.classification.each do |taxon|
      puts classification[i] + "\t" + taxon
      i += 1
    end
    puts classification[i] + "\t" + @gbkObj.organism
  end


end                             # End of Class GenbankManipulator


