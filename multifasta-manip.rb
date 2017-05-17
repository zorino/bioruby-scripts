#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# review:  	
# date:    	12-11-25
# version: 	0.01
# licence:  	

require_relative 'multifasta-manip-class.rb'


# Function that 

# Fct: Return unique sequence
def uniqueSequences file

  seqIDs = {}
  printSeq = 1
  File.open(file,"r") do |f|
    while l = f.gets
      if l[0] == ">"
        # key = l.gsub("[","").gsub("]","")
        key = l.split("\n")[0].split(" ")[0]
        if seqIDs.has_key? key
          printSeq = 0
        else
          seqIDs["#{key}"] = 0
          printSeq = 1
          puts l
        end
      elsif printSeq == 1
        puts l
      end
    end
  end

end



# Main #

@usage = 
  "
  multifasta-manip.rb [options] <multifasta>

	options :	geneNames
			getSeq <regex>
			split (will split each sequence)
			splitNb <number of file output> 
			merge <number of N's>
			length <length>
			sortLength <asc||dsc> <out>
			getseq-loc <seqName> <beg..end> <strand>
			getseq-unique
"

abort "No input file ! \n #{@usage}" if ARGV[0].nil?

fasta = FastaParser.new(ARGV[ARGV.length-1])
opt = ARGV[0]

case opt.downcase

when 'genenames'
  fasta.printGenes

when 'getseq'
  abort "You didn't enter a Gene Name or ID" if ARGV.length <= 2
  name = ARGV[1]
  gene = fasta.getGene(name)
  puts gene

when 'split'
  fasta.split
  
when 'splitnb'
  fasta.splitNumber ARGV[1].to_i

when 'merge'
  fasta.merge ARGV[1].to_i

when 'length'
  fasta.lengthExtract ARGV[1]

when 'sortlength'
  fasta.sortLength ARGV[1], ARGV[2]

when 'getseq-loc'
  if ARGV.length < 4
    abort "You need to provide sequence name, beginning, end and strand"
  end
  fasta.getSeqLoc ARGV[1], ARGV[2], ARGV[3]

when 'getseq-unique'
  uniqueSequences ARGV[1]
else
  puts @usage

end
