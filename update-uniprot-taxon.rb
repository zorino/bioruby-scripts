#!/usr/bin/env ruby
# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# review:  	
# date:    	13-05-08
# version: 	1.0
# licence:  	

require 'getopt/long'
require 'net/http'
require 'uri'

usage = "

usage : update-uniprot-taxon.rb --taxon(-t) <Taxon Root Number>

	Options : 	--new (for new releases)
			--leaf <files of taxon> (only download those taxa)

"

if ARGV.length < 1
  puts "\nWill execute the update with default parameters : \n\n \t Root Taxon = cellular organism and Not like a New version"
  puts usage
end


# OPTIONS Settings
include Getopt

begin
  opt = Long.getopts( ["--taxon","-t", Getopt::REQUIRED],
                      ["--new", "-n", Getopt::BOOLEAN],
                      ["--leaf", "l", Getopt::REQUIRED])
rescue
  puts usage
  exit 0
end

new = 0
if opt["taxon"]
  taxon = opt["taxon"]
  if opt["new"]
    new = 1
  end
else
  taxon = "131567" # cellular organism
end

if opt["leaf"]
  taxonH = {}
  File.open(opt["leaf"]) do |f|
    while l = f.gets
      l.chomp!
      taxonH[l.to_s] = "0"
    end
  end
end
# End of options setting
# The rest of the main is at the end of the file


#-# Function #-#

# Fetch taxons from uniprot
def getResponse url
  code = ""
  body = ""
  while code != "200" or body.include? "DOCTYPE HTML" 
    begin
      uriList = URI.parse(url)
      http = Net::HTTP.new(uriList.host,uriList.port)
      http.read_timeout = 500
      request = Net::HTTP::Get.new(uriList.request_uri)
      response = http.request(request)
      code = response.code.to_s
      body = response.body.to_s
    rescue
      puts "Problem with the Connection !"
      puts "Will wait 25 before retry .. "
      sleep 25
      retry
    end
 end
  return body
end

# Format sequences as fasta (60 widht)
def string2seq seq
  numN = (seq.length-1) / 60
  for i in 1..numN
    index = i*61 - 1
    seq.insert(index,"\n")
  end
  return seq 
end


# Add sequences to new sequences to treat
def addNewSeq header, seq
  time = Time.new
  timeS = "#{time.year}_#{time.month}"
  newSeq = string2seq seq  
  File.open("000-#{timeS}-New-Sequences.fasta",'a') do |newF|
    newF << "#{header}\n"
    newF << "#{newSeq}\n"
  end
end

 
# Check if taxon exist
def taxonExists? oldtaxon, seq, seqName
  if oldtaxon.has_key?("#{seqName}")
    if seq.to_s != oldtaxon["#{seqName}"].to_s
      # puts "---------------= VS =------------------"
      addNewSeq(seqName,seq)
    end
  else
    # puts "oldtaxon has not that key"
    addNewSeq(seqName,seq)
  end
end


#-# Main #-# 

taxonsList = getResponse("http://www.uniprot.org/taxonomy/?query=ancestor:#{taxon}&format=list")
taxons = taxonsList.to_s.split("\n")
time = Time.new
timeS = "#{time.year}_#{time.month}"

taxons.each do |id|

  if opt["leaf"] and ! taxonH.has_key? id.to_s
    next     # Not a leaf taxa
  end

  url2Fetch = "http://www.uniprot.org/uniprot/?query=taxonomy%3a#{id}&force=yes&format=fasta"
  puts "check #{id}"

  if File.exist?("#{id}.fasta") and new == 1 # File (taxon) exists, checking for NEW seq

    puts "fetch #{id}"
    seqs = getResponse(url2Fetch)

    if seqs != ""

      oldtaxon = {}
      
      # Open old file and fill hash of prots
      File.open("#{id}.fasta", "r") do |oldf|
        protName = ""
        seq = ""
        while l = oldf.gets
          l.chomp!
          if l[0] == ">"
            oldtaxon["#{protName}"] = seq unless protName == ""
            protName = l.to_s
            seq = ""
          else
            seq = "#{seq}" + l.to_s
          end
        end
        oldtaxon["#{protName}"] = seq unless protName == ""
      end
      
      puts "Check for new sequences in file"
      # Check for new sequences in file
      seqArray = seqs.split("\n")
      seq = ""
      seqName = ""
      for i in 0..seqArray.length
        if i == seqArray.length
          taxonExists?(oldtaxon, seq, seqName)
          # puts seq
        elsif  seqArray[i][0] == ">"
          taxonExists?(oldtaxon, seq, seqName) unless seqName == ""
          seqName = seqArray[i].to_s.gsub("\n","")
          seq = ""
        else
          seq = "#{seq}" + seqArray[i].to_s 
        end
      end
      
    end
    
  elsif  File.exist?("#{id}.fasta") and new == 0 # File (taxon) exists, not checking for NEW seq
    next
  else                    	# File doesn't exist .. a new taxon (& dump all in NEW)
    puts "fetch #{id}"
    seqs = getResponse(url2Fetch)
    if seqs != ""
      File.open("000-#{timeS}-New-Sequences.fasta","a") do |newF|
        newF.write(seqs)
      end
    end
  end                           # End If file exists
  
  # Save New File Taxon File and Overwrite last one
  File.open("#{id}.fasta", "w") {|f|  f.write(seqs.to_s) } unless seqs == ""
  
end	# End of #+# Main #+#
