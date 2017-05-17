#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# review:  	
# date:    	14-09-16
# version: 	0.01
# licence:  	


require 'bio'


if ARGV.length < 3
  abort "usage: tRNAscan-to-gbk.rb <gbk file> <tRNAscan output> <prefix>"
end

genbankFile = ARGV[0]
tRNAFile = ARGV[1]
prefix = ARGV[2]


def loadRNAfile tRNAFile

  inResults = 0
  results = []

  File.open(tRNAFile,"r") do |f|
    while l = f.gets

      if inResults != 0
        lA = l.split("\t")
        results.push( {
          tRNANb: lA[1],
          seqBeg: lA[2].strip,
          seqEnd: lA[3].strip,
          tRNAtype: lA[4],
          anticodon: lA[5],
          covScore: lA[6]
        })
      end

      if l[0..7] == "--------"
        inResults = 1
      end

    end    
  end

  results.sort_by{ |h| h[:seqBeg].to_i }

end 


def addFeatureToGenbank genbankFile, tRNArray, prefix

  gbk = Bio::GenBank.new(File.read(genbankFile))

  index = 1

  tRNArray.each do |rna|

    if rna[:seqBeg] > rna[:seqEnd]
      feature = Bio::Feature.new('tRNA',"complement(#{rna[:seqEnd]}..#{rna[:seqBeg]})")
    else
      feature = Bio::Feature.new('tRNA',"#{rna[:seqBeg]}..#{rna[:seqEnd]}")
    end

    tRNANb = format('%02d', index)
    feature.append(Bio::Feature::Qualifier.new('locus_tag', "#{prefix}_tRNA#{tRNANb}"))
    feature.append(Bio::Feature::Qualifier.new('product', "tRNA-#{rna[:tRNAtype]}"))

    gbk.features.push(feature)

    index += 1

  end

  puts gbk.to_biosequence.output(:genbank)

end

tRNArray = loadRNAfile tRNAFile
addFeatureToGenbank genbankFile, tRNArray, prefix
