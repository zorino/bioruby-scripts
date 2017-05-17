#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# review:  	
# date:    	12-12-14
# version: 	0.01
# licence:  	

# MSA parser between format

require 'bio'


usage = 
"
   clustalMSAparser.rb <clustal aln> <out format> <out file>  
		outfmt values 
				clustal for ClustalW
				fasta for FASTA
				phylip 
				phylipnon
				msf for MSF
				molphy for Molphy
"

fclu = ARGV[0] or abort usage
outfmt = ARGV[1] or abort usage
fout = ARGV[2] or abort usage

# Reads in a ClustalW-formatted multiple sequence alignment
# from a file named "infile_clustalw.aln" and stores it in 'report'.
report = Bio::ClustalW::Report.new(File.read(fclu))

# Accesses the actual alignment.
msa = report.alignment


# Write the output

File.open("#{fout}", "w") do |f|
  f.write(msa.output(:"#{outfmt}"))
end
