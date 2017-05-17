#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# review:  	
# date:    	14-02-23
# version: 	0.01
# licence:  	


# Define the usage of this script
usage = "\n$ cd-hit-parser.rb <cdHitFile> <samplesFile>"
usage += "\n\n  ! Results header should look like this to map the samples file"
usage += "\n   >FM209186_PLES_28951"
usage += "\n   >SeqID_LocusTag\n\n"

if ARGV.length < 2
  abort usage
end

cdHitFile = ARGV[0]
samplesFile = ARGV[1]


# Load all sample in a hash first
samplesHash = {}
header = ""
File.open(samplesFile,"r") do |f|
  header = "Clusters"
  i = 1
  while l=f.gets
    header += "\t#{l.chomp!}"
    samplesHash[l.chomp] = i
    i += 1
  end
  header += "\n"
end

puts header

lineArray = Array.new((samplesHash.length+1),"-")
# Iterate of cd-hit files and print in a matrix tsv file
File.open(cdHitFile,"r") do |f|
  while l=f.gets
    if l[0] == ">"              # cluster id
      if lineArray[0] != "-"
        puts lineArray.join("\t")
      end
      lineArray.fill("-")
      clusterId =  l.chomp![1..-1]
      lineArray[0] = clusterId
    else                        # member of cluster
      lA = l.chomp!.split(", >")
      aaNumber = lA[0].split()[1].gsub("aa,","").to_i
      hit = lA[1].split("... ")
      sampleId = hit[0].split("_")[0]
      locus = hit[0].split("_")[1..-1].join("_")
      idPercentage = ""
      if ! hit[1].match(/\*/).nil?
        idPercentage = "100"
      else
        idPercentage = hit[1].gsub("at ","").gsub("%","")
      end
      # puts "#{aaNumber}\t#{sampleId}\t#{locus}\t#{idPercentage}"
      # puts sampleId
      index = samplesHash[sampleId]
      # puts index
      lineArray[index] = "#{locus} \[#{idPercentage}\]"
    end
  end
end
