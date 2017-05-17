#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# review:  	
# date:    	14-01-13
# version: 	0.01
# licence:  	

require 'bio'

if ARGV.length < 2
  puts " usage : prodigal-add-annotation.rb <annotation sum file> <locus product name> <genbank files dir>"
  abort
end

annotationFile = ARGV[0]
nameAssociation = ARGV[1]
gbkFiles = ARGV[2]

annotations = {}

File.open(annotationFile,"r") do |f|
  while l = f.gets
    lA = l.chomp!.split("\t")    
    if ! annotations.has_key? lA[0]
      annotations[lA[0]] = lA[1]
    end
  end
end

productName = {}
File.open(nameAssociation,"r") do |f|
  while l = f.gets
    lA = l.chomp!.split("\t")    
    if ! productName.has_key? lA[0]
      productName[lA[0]] = lA[1]
    end
  end
end


# # Main LOOP for each gbk File # #

Dir.glob(gbkFiles).each do |gbkFile|

  puts gbkFile
  contig = gbkFile.split("/")[-1].split(".")[0]

  puts contig

  gbk = Bio::GenBank.new(File.open(gbkFile,"r").read)

  i = 0
  gbk.each_cds do |cds|

    ftArray = cds.qualifiers

    i += 1
    hit = annotations[contig+"_"+i.to_s]

    if hit != nil

      if hit.length > 90
        hit.insert(88," ")
      end

      hitArray = hit.split("|")

      locus = hitArray[3]
      gene = hitArray[4]
      product = productName[locus]

      if gene != ""
        qGene = Bio::Feature::Qualifier.new('gene', gene)
        ftArray.push(qGene)
      end

      qProd = Bio::Feature::Qualifier.new('product', product)
      qNote = Bio::Feature::Qualifier.new('note', "correspond to #{locus} locus")
      ftArray.push(qProd,qNote)

      cds.qualifiers = ftArray

    end

  end

  File.open("#{contig}-new.gbk", "w") do |f| 
    f.write(gbk.to_biosequence.output(:genbank))
  end

  # "NC_002516|complement(5919585..5920715)|NP_253945.1|PA5258||hypothetical protein "

end

