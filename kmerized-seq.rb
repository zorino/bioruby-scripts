#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# review:  	
# date:    	13-12-16
# version: 	0.01
# licence:  	

require 'bio'
require 'json'
require 'google_hash'

def getGenbankAnnotationForKmer kmerPosition, gbk

  position = []
  features = []

  gbk.features.each do |ft|

    looping = 1

    next if ft.feature.to_s != "CDS" and ! ft.feature.to_s.include? "RNA"

    positionS = ft.position.gsub("complement","").gsub("(","").gsub(")","")

    position = positionS.split("..")

    if kmerPosition.to_i >= position[0].to_i and kmerPosition.to_i <= position[1].to_i
    
      kmerFeature = {}

      gene = ""
      product = ""
      locus_tag = ""
      protein_id = ""

      ft.qualifiers.each do |q|

        case q.qualifier
        when "gene"
          gene = q.value
        when "product"
          product = q.value
        when "locus_tag"
          locus_tag = q.value
        when "protein_id"
          protein_id = q.value
        end
 
      end

      kmerFeature[:gene] = gene
      kmerFeature[:product] = product
      kmerFeature[:locus_tag] = locus_tag
      kmerFeature[:protein_id] = protein_id
      kmerFeature[:posBeg] = position[0].to_i
      kmerFeature[:posEnd] = position[1].to_i

      features.push(kmerFeature)

    elsif kmerPosition.to_i < position[0].to_i

      looping = 0

    end

    break if looping == 0

  end

  return features

  # if features.empty? 
  #   return {ft: features, position: {posBeg: '0', posEnd: '0'}}
  # else
  #   return {ft: features, position: {posBeg: position[0], posEnd: position[1]}}
  # end

end



def kmerizeFasta(kmersHash,sequence,kmerSize,gbk)

  i = 0
  subSeq = ""
  # kmers = []
  # position = {posBeg: '-1', posEnd: '-1'}

  kmerFeatures = [{posBeg: 0, posEnd: 0}]

  # k-mer nucleotide sequence; coverage value; first nucleotide of parents; last nucleotide of children
  while i < (sequence.length-kmerSize.to_i+1)

    if kmerFeatures.empty?
      kmerFeatures = getGenbankAnnotationForKmer i+1, gbk
    elsif ! (i+1.to_i >= kmerFeatures[0][:posBeg].to_i and i+1.to_i <= kmerFeatures[0][:posEnd].to_i)
      kmerFeatures = getGenbankAnnotationForKmer i+1, gbk   
    end

    puts "#{i} | #{kmerFeatures}"

    subSeq = sequence[i..(i+kmerSize.to_i-1)]
    parentArray = []
    childArray = []
    parent = ""
    child = ""
    parent = sequence[i-1] if i>0
    child = sequence[i+kmerSize.to_i] if (i+kmerSize.to_i) < sequence.length

    # Check for kmer key in hash .. might be in a second pass

    # if kmersHash.has_key? subSeq
    #   parentArray = kmersHash[subSeq][:parent]
    #   childArray = kmersHash[subSeq][:child]
    #   parentArray.push(parent) if parent != ""
    #   childArray.push(child) if child != ""
    #   depth = kmersHash[subSeq][:depth] + 1
    #   kmersHash[subSeq] = {depth: depth, parent: parentArray, child: childArray, features: kmerFeatures}
    # else
      parentArray.push(parent) if parent != ""
      childArray.push(child) if child != ""
      kmersHash[subSeq] = {depth: 1, parent: parentArray, child: childArray, features: kmerFeatures}
    # end

    # kmers.push("#{subSeq};10;#{parent};#{child}") 1
    i += 1
  end

  # return kmers

end



# # MAIN # # 

usage = "
	kmerized-seq.sh <input-fasta> <kmer size> <genbank file>
"

if ARGV.length < 2
  abort usage
else
  seqFile = ARGV[0]
  kmerSize = ARGV[1]
  flatfile = Bio::FlatFile.auto(seqFile)
  if flatfile.dbclass.to_s != "Bio::FastaFormat"
    abort "#{seqFile} is not a Fasta File !!"
  end
  if ARGV[2]
    gbk = Bio::GenBank.new(File.read(ARGV[2]))    
  end
end

# fOut = File.open("kmers.txt","w")

nbOfSeq = 0
nbOfKmers = 0
# kmersHash = {}
kmersHash = GoogleHashDenseRubyToRuby.new

flatfile.each_entry do |fasta|
  seq = fasta.seq
  # kmers = kmerizeFasta(seq.upcase,kmerSize)
  kmerizeFasta(kmersHash,seq.upcase,kmerSize,gbk)
  # fOut.write(kmers.join("\n")+"\n")
  nbOfKmers += kmersHash.length
  nbOfSeq += 1
end

File.open("kmers.txt","w") do |f|
  kmersHash.each_key do |k|
    d = kmersHash[k][:depth]
    p = kmersHash[k][:parent].join(" ")
    c = kmersHash[k][:child].join(" ")
    ft = kmersHash[k][:features].join(",") 
    f << "#{k};#{d};#{p};#{c};#{ft}\n"
  end
end



# fOut.close

puts "#{nbOfSeq} kmerized sequences\n"
puts "#{nbOfKmers} total kmers of #{kmerSize} length\n"

