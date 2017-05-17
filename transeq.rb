#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# review:  	
# date:    	12-12-15
# version: 	0.01
# licence:  	

require 'bio'

usage =
"
   transeq.rb <in.fa> [strand]
"

strand = 1
if ARGV[1] == "0" or ARGV[1] == "-1"
  strand = 0
end

ntF = ARGV[0] or abort usage

ntSeq = Bio::FlatFile.auto(ntF)

ntSeq.each_entry do |e|
  # pep = Bio::Sequence.new(e.naseq.reverse_complement!.translate)
  pep = Bio::Sequence.new(e.naseq.translate)
  if strand == 0
    pep = Bio::Sequence.new(e.naseq.reverse_complement!.translate)
  end
  puts pep.output_fasta(e.definition,60)
end
