#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# review:  	
# date:    	14-09-17
# version: 	0.01
# licence:  	

require 'bio'


def loadBlast8File blast8File

  # Query should be what you want to identified in the sample
  # AP013063|210161..211701||SM39_rRNA01||16S       SM909   99.74   1541    4       0       1       1541    301164  302704  0       1945.9

  results = []

  File.open(blast8File,"r") do |f|

    while l = f.gets

      lA = l.split("\t")

      hit = lA[0].split("|")
      location = hit[1].gsub("complement","").gsub(")","").gsub("(","").split("..")
      hitLocusTag = hit[3]
      hitProduct = hit[5]
      
      pId = lA[2]
      alignmentLength = lA[3]
      seqBeg = lA[8]
      seqEnd = lA[9]

      hitLength = location[1].to_i - location[0].to_i + 1

      if alignmentLength.to_i == hitLength and pId.to_f > 90
 
        to_add = 1

        results.each do |res|

          if res[:seqBeg] == seqBeg and res[:seqEnd] == seqEnd
            to_add = 0
            if res[:pId].to_f >= pId.to_f
              res[:pId] = pId.to_f
              res[:hitLocusTag] = hitLocusTag
              res[:hitProduct] = hitProduct
              break
            end
          end

        end

        if to_add == 1
          results.push({ hitLocusTag: hitLocusTag,
                         hitProduct: hitProduct,
                         pId: pId,
                         seqBeg: seqBeg,
                         seqEnd: seqEnd
                       })
        end

      end

    end    
  end

  results.sort_by{ |h| h[:seqBeg].to_i }

end


def addFeatureToGenbank genbankFile, rRNArray, prefix

  gbk = Bio::GenBank.new(File.read(genbankFile))

  index = 1

  rRNArray.each do |rna|

    if rna[:seqBeg] > rna[:seqEnd]
      feature = Bio::Feature.new('rRNA',"complement(#{rna[:seqEnd]}..#{rna[:seqBeg]})")
    else
      feature = Bio::Feature.new('rRNA',"#{rna[:seqBeg]}..#{rna[:seqEnd]}")
    end
    


    tRNANb = format('%02d', index)

    feature.append(Bio::Feature::Qualifier.new('locus_tag', "#{prefix}_rRNA#{tRNANb}"))
    feature.append(Bio::Feature::Qualifier.new('product', "#{rna[:hitProduct]} ribosomal RNA"))

    gbk.features.push(feature)

    index += 1

  end

  puts gbk.to_biosequence.output(:genbank)

end


# # MAIN # #

genbankFile = ARGV[0]
blast8File = ARGV[1]            # purge the results before
prefix = ARGV[2]

results = loadBlast8File blast8File
addFeatureToGenbank genbankFile, results, prefix


