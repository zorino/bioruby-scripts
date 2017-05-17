#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# review:  	
# date:    	13-06-06
# version: 	0.01
# licence:  	

require 'bio'

usage = "Usage : $ remoteblast.rb <fasta> [blastType]"


if ARGV[0]
  seqFile = Bio::FlatFile.auto(ARGV[0])
else
  abort usage
end


@blastType = "blastp"

if ARGV[1]
  @blastType = ARGV[1].to_s
end


def runBlast sequence_text
  # To run an actual BLAST analysis:
  #   1. create a BLAST factory
  blast_factory = Bio::Blast.remote(@blastType, 'nr-aa',
                                    '-e 0.0001', 'genomenet')

  resultArray = Array.new

  #   2. run the actual BLAST by querying the factory
  report = blast_factory.query(sequence_text)

  report.each_hit do |hit|
    resultArray.push("#{hit.query_start}\t#{hit.query_end}\t#{hit.definition}\t#{hit.evalue}")
    # puts hit.definition
    # puts hit.evalue
    # puts ""
    break
  end

  resultArray

end

if seqFile.dbclass.to_s != "Bio::FastaFormat"
  abort "Wrong File Format .. need a fasta file"
end


seqFile.each_entry do |s|
  runBlast(s.seq).each do |res|
    puts s.definition.chomp + "\t#{res}"
  end
end


