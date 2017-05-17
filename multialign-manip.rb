#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# review:  	
# date:    	14-03-01
# version: 	0.01
# licence:  	

require 'bio'
require 'bio-alignment'
require 'bio-alignment/bioruby'

include Bio::BioAlignment

usage = " 
     multialign-manip.rb  <options>

	opt:  	concatenate   	<multialign files> # Need to have the same genes order & length in each MSA file
		addGapSeq	<samplesOrder txt> <seqLength> <multalign file>
		consensus	<multialign file>
"

def addGapSeq seqFile, seqLen, multialignFile

  seqArray = []
  File.open(seqFile,"r") do |f|
    while l=f.gets
      seqArray.push(l.chomp!)
    end
  end

  File.open(multialignFile,"r") do |f| 

    i = 0
    nextSeq = ""

    while l = f.gets
      if l[0] == ">"
        currentSeq = l.chomp!.gsub(">","")
        if currentSeq == seqArray[i]
          puts l  
        else
          while seqArray[i] != currentSeq
            bioSeq = Bio::Sequence.auto("-"*seqLen.to_i)
            puts bioSeq.output_fasta("#{seqArray[i]}", 60)
            i+=1
          end
          puts l
        end
        i += 1
      else
        puts l
      end
    end

    while i < (seqArray.length - 1)
      bioSeq = Bio::Sequence.auto("-"*seqLen.to_i)
      puts bioSeq.output_fasta("#{seqArray[i]}", 60)
      i += 1
    end

  end

end



def concatenate geneList, multiFastaFile

  def addSeq flat, samplesHash

    flat.each_entry do |ent|
      if samplesHash.has_key? ent.entry_id.to_s
        samplesHash[ent.entry_id] += ent.seq.to_s
      else
        samplesHash[ent.entry_id] = ent.seq.to_s
      end
    end

  end

  samplesArray = []

  File.open(geneList,"r") do |f|
    while l = f.gets      
      if l.chomp != ""
        samplesArray.push(l.chomp!)
      end
    end
  end

  samplesHash = {}

  multiFastaFile.each do |f|
    flat = Bio::FlatFile.auto(f)
    addSeq flat, samplesHash
  end


  samplesHash.each do |k,v|
    outSeq = Bio::Sequence.new(v)
    puts outSeq.output_fasta(k)
  end

  # for i in 0..samplesArray.length-1
  #   seqString = ""    

  #   multiFastaFile.each do |f|
  #     flat = Bio::FlatFile.auto(f)
  #     seqString += addSeq flat, samplesArray[i].to_s
  #   end

  #   outSeq = Bio::Sequence.new(seqString)
  #   puts outSeq.output_fasta(samplesArray[i].to_s)

  # end


end



def consensus alignmentFile


  aln = Alignment.new
  # puts aln.methods

  Bio::FlatFile.auto(alignmentFile).each_entry do |entry|
    aln << entry
  end

  bioruby_aln = aln.to_bioruby_alignment

  # puts bioruby_aln.methods

  puts bioruby_aln.consensus_string

  bioruby_aln.each_site do |site|
    
    
  end


end


# # Main # #


if ARGV.length < 1
  abort usage
end

case ARGV[0].downcase

when 'concatenate'
  concatenate ARGV[1], ARGV[2..-1]
when 'addgapseq'
  addGapSeq ARGV[1], ARGV[2], ARGV[3]
when 'consensus'
  consensus ARGV[1]
else
  puts usage

end

