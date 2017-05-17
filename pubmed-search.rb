#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# review:  	
# date:    	12-12-04
# version: 	0.01
# licence:  	

require 'bio'

if ARGV.length < 1
  abort "pubmed-search.rb 'keywords' [medline|bibtex|cite]"
end

# search and print result
Bio::NCBI.default_email = 'default@default.com'

entries = Bio::PubMed.esearch(ARGV[0], { "retmax" => 100000 } )

Bio::PubMed.efetch(entries).each do |entry|
  medline = Bio::MEDLINE.new(entry)
  reference = medline.reference

  if ARGV[1] == "medline"
    puts entry
  elsif ARGV[1] == "bibtex"
    puts reference.bibtex
  elsif (ARGV[1] == "cite" or ARGV[1] == "citation")
    puts reference.format
  else
    puts entry
  end

  puts "//"
end


# Bio::PubMed.search("(genome AND analysis) OR bioinformatics").each do |x|
#   p x
# end

# # To retrieve the MEDLINE entry for a given PubMed ID:
# puts Bio::PubMed.efetch("10592173", "14693808")
# puts Bio::PubMed.query("10592173")
# puts Bio::PubMed.pmfetch("10592173")

# # This can be converted into a Bio::MEDLINE object:
# manuscript = Bio::PubMed.query("10592173")
# medline = Bio::MEDLINE.new(manuscript)
